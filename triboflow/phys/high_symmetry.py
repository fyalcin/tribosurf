#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:52:15 2020

Calculate the High Simmetry (HS) points for slab and interface

@author: gl
"""

import numpy as np
from monty.json import jsanitize
from ase import Atoms
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.structure_matcher import StructureMatcher
from triboflow.utils.structure_manipulation import StackAlignedSlabs, \
    CleanUpSiteProperties, ReCenterAlignedSlabs

# =============================================================================
# CALCULATE THE HS POINTS FOR A SLAB
# =============================================================================


def GetSlabHS(slab, allowed_sites=['ontop', 'bridge', 'hollow'], to_array=False): 
    """
    Calculate the High Simmetry (HS) points for a material provided as input,
    which need to be a pymatgen object (Slab or Structure).
    It returns the unique HS sites for the selected slab and all the sites
    unfolded inside the lattice cell of the object.

    Parameters
    ----------
    slab : pymatgen.core.surface.Slab (pymatgen.core.structure.Structure)
        The slab of which you want to calculate the surface HS points.
        The input need to be a slab, i.e. an atomic structure non-periodic
        along the z-direction. However, type(slab) could be either a Slab
        object or a Structure object.
        
    allowed_sites : list, optional
        The type of the HS sites that you want to extract from your slab.
        Unless specific reasons, just leave the default and work on the sites 
        of your interest later. The default is ['ontop', 'bridge', 'hollow'].
        
    to_array : bool, optional
        If you want to return arrays or lists. The default is False.

    Returns
    -------
    hs : dict
        Contain the unique surfacial HS sites of slab. Multiple sites are
        enumerated to distinguish them. To be coherent with the notation, 
        even single sites of one type are renamed in the same way.
        Example:
        - Slab has one 'bridge' site  -> named as 'bridge_1'.
        - Slab has two 'bridge' sites -> named as 'bridge_1', 'bridge_2'.
        
    hs_all : dict 
        Surfacial HS sites that has been unfolded (replicated) across the whole
        lattice cell of slab. Sites are renamed in the same way as hs.
        It is normally easy to obtain the HS sites for a replicated defect-free
        unit cell, but this operation may become tricky for a slab that has 
        been strained, rotated and cut to be matched with another one to form 
        an hetero interface.

    """

    adsf = AdsorbateSiteFinder(slab)
    
    # Extract the unique HS points for the given surface in the unit cell
    unique_sites = adsf.find_adsorption_sites( distance = 0, 
                                               symm_reduce = 0.01, 
                                               near_reduce = 0.01, 
                                               no_obtuse_hollow = True )
    
    #Extract all the HS points for the given surface in the lattice cell
    replica_sites = adsf.find_adsorption_sites( distance = 0, 
                                                symm_reduce = 0, 
                                                near_reduce = 0.01, 
                                                no_obtuse_hollow = True )
    
    # Identify the unique HS points of the slab and rename them
    hs = {} 
    for key in allowed_sites:
        if unique_sites[key] != []:
            for n, data in enumerate(unique_sites[key]):
                hs[key+'_'+str(n+1)] = data
            
    # Recognize of which unique HS point each site is the replica
    hs_all = {}
    for k in hs.keys():
        hs_all[k] = []   
    for key in hs_all.keys():
        for site in replica_sites['all']:         
            pts_to_evaluate = [hs[key].copy()]
            pts_to_evaluate.append(np.array(site))
            pts = adsf.symm_reduce(pts_to_evaluate)     
            if len(pts) == 1:
                hs_all[key].append(site)
    
    # Add a key to the dictionary related to all the sites. Not necessary
    #hs['all'] = unique_sites['all']
    #hs_all['all'] = replica_sites['all']
    
    # Convert the elements of the dictionaries to proper numpy arrays and 
    # remove the z coordinate of the HS sites
    hs = NormalizeHSDict(hs, to_array)
    hs_all = NormalizeHSDict(hs_all, to_array)

    return hs, hs_all


def HS_DictConverter(hs, to_array=True):
    """
    Modify the type of the elements of the HS dictionary to list or np.ndarray.

    Parameters
    ----------
    hs : dict
        Dictionary containing the High Symmetry points.
    to_array : bool, optional
        If set to True convert to array, otherwise convert to list. 
        The default is True.

    Raises
    ------
    ValueError
        Raised if the dictionary values are of different types.
        Print to stdout: "Your dictionary is weird, values have mixed types"

    Returns
    -------
    hs_new : dict
        New HS dictionary converted to the desired type.

    """
    
    hs_new = {}
    dict_types = list( set(type(k) for k in hs.values()) )
    
    try: 
        assert(len(dict_types) == 1)
        
        typ = dict_types[0]
        if to_array:
            if typ == list:
                for k in hs.keys():
                    hs_new[k] = np.array(hs[k])    
            else:
                return hs
            
        else:
            if typ == np.ndarray:
                for k in hs.keys():
                    hs_new[k] = hs[k].tolist() 
            else:
                return hs
            
        return hs_new
            
    except:
        raise ValueError('Your dictionary is weird, values have mixed types')


# =============================================================================
# CALCULATE HS POINTS FOR AN INTERFACE
# =============================================================================


def GetInterfaceHS(hs_1, hs_2, cell, to_array=False, z_red=True):
    """
    Calculate the HS sites for a hetero interface by combining the HS sites of
    the bottom slab (hs_1) with the upper slab (hs_2) 

    Parameters
    ----------
    hs_1 : dict
        High Symmetry points of the bottom slab
    hs_2 : dict
        High Symmetry points of the upper slab
    to_array : bool, optional
        If set to True return an HS dictionary containing arrays, else lists.
        The default is False.
    z_red : bool, optional
        Remove the z-coordinates from the translations. The default is True.

    Returns
    -------
    hs : dict
        High Symmetry points of the hetero interface.

    """
    
    hs = {}
    
    typ_1 = list( set(type(k) for k in hs_1.values()) )[0]
    if typ_1 == list:
        hs_1 = HS_DictConverter(hs_1, to_array=True)
        
    typ_2 = list( set(type(k) for k in hs_1.values()) )[0]
    if typ_2 == list:
        hs_2 = HS_DictConverter(hs_2, to_array=True)
        
    
    # Calculate the shift between each HS point of the first material with each
    # HS point of the second material
    for k1, v1 in hs_1.items():
        for k2, v2 in hs_2.items():  
            shifts_stack = []
            for el_d1 in v1:
                shifts_stack.append( v2 - el_d1 )
                
            hs[k1+'-'+k2] = np.concatenate(shifts_stack, axis=0)
    
    hs = PBC_HSPoints(hs, cell, to_array=to_array, z_red=z_red)
        
    return hs


# =============================================================================
# TOOLS FOR HS DICTIONARIES
# =============================================================================


def FixHSDicts(hs_unique, hs_all, top_aligned, bot_aligned,
               ltol=0.01, stol=0.01, angle_tol=0.01,
               primitive_cell=False, scale=False):
    """Remove duplicate shifts from the hs points and assigne the replicas correctly.
    
    A StructureMatcher is defined with the selected tolerances and options and
    then used to remove equivalent shifts from the high-symmetry points
    dictionaries and assign the replicated points correctly to their unique
    counterparts using the <RemoveEquivalentShifts> and <AssignReplicatePoints>
    functions.
    

    Parameters
    ----------
    hs_unique : dict
        Unique high-symmetry points of the interface from <GetInterfaceHS>.
    hs_all : dict
        Replicated high-symmetry points of the interface from <GetInterfaceHS>.
    top_aligned : pymatgen.core.surface.Slab or pymatgen.core.structure.Structure
        The top slab of the interface
    bot_aligned : pymatgen.core.surface.Slab or pymatgen.core.structure.Structure
        The bottom slab of the interfaces
    ltol : float, optional
       Fractional length tolerance. The default is 0.01.
    stol : float, optional
        Site tolerance. The default is 0.01.
    angle_tol : float, optional
        Angle tolerance in degrees. The default is 0.01.
    primitive_cell : bool, optional
        If true: input structures will be reduced to primitive cells prior to
        matching. The default is False.
    scale : bool, optional
        Input structures are scaled to equivalent volume if true; For exact
        matching, set to False. The default is False.

    Returns
    -------
    c_u : TYPE
        DESCRIPTION.
    c_all : TYPE
        DESCRIPTION.

    """
    top_slab, bot_slab = ReCenterAlignedSlabs(top_aligned,
                                              bot_aligned,
                                              d=4.5)
    
    struct_match = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol,
                                primitive_cell=primitive_cell, scale=scale)
    
    #Use the structure matcher to find shifts leading to equivalent interfaces
    #and pop these entries out of the dictionaries.
    c_u, c_a = RemoveEquivalentShifts(hs_u=hs_unique.copy(),
                                      hs_a=hs_all.copy(),
                                      top_slab=top_slab,
                                      bot_slab=bot_slab,
                                      structure_matcher=struct_match)
    
    c_all = AssignReplicatePoints(hs_u=c_u,
                                   hs_a=c_a,
                                   top_slab=top_slab,
                                   bot_slab=bot_slab,
                                   structure_matcher=struct_match)
    
    return c_u, c_all
    

def AssignReplicatePoints(hs_u, hs_a, top_slab, bot_slab, structure_matcher):
    """Assign the replicated high-symmetry points to the correct unique ones.
    
    Although most of the high-symmetry points should be assigned to the correct
    lable, there is the occasional shift that is equivalent for two lables.
    This function imploys the StructureMatcher to match the replicated points
    to their unique counterparts, so the energy can later be transfered
    correctly.
    

    Parameters
    ----------
    hs_u : dict
        Unique high-symmetry points of the interface.
    hs_a : dict
        All high-symmetry points of the interface.
    top_slab : pymatgen.core.surface.Slab or pymatgen.core.structure.Structure
        The top slab of the interface
    bot_slab : pymatgen.core.surface.Slab or pymatgen.core.structure.Structure
        The bottom slab of the interfaces
    structure_matcher : pymatgen.analysis.structure_matcher.StructureMatcher
        Class to find equivalent structures (mirrors, rotations, etc...)

    Returns
    -------
    new_hsp_dict_a : dict
        All high Symmetry points of the interface without duplicated entries.

    """
    all_shifts = []
    for key, value in hs_a.items():
        if all_shifts == []:
            all_shifts = value
        else:
            all_shifts = np.concatenate([all_shifts, value], axis=0).tolist()
    all_shifts = np.unique(all_shifts, axis=0)

    
    new_hsp_dict_a = {}
    for key, value in hs_u.items():
        unique_struct = StackAlignedSlabs(bot_slab,
                                          top_slab,
                                          top_shift = [value[0][0],
                                                       value[0][1],
                                                       0])
        unique_struct = CleanUpSiteProperties(unique_struct)
        for shift in all_shifts:
            test_struct = StackAlignedSlabs(bot_slab,
                                            top_slab,
                                            top_shift = [shift[0],
                                                         shift[1],
                                                         0])
            test_struct = CleanUpSiteProperties(test_struct)
            if structure_matcher.fit(unique_struct, test_struct):
                new_hsp_dict_a.setdefault(key, []).append(shift)
    
    return new_hsp_dict_a

def RemoveEquivalentShifts(hs_u, hs_a, top_slab, bot_slab, structure_matcher):
    """Remove equivalent shifts from an interface high-symmetry point dict.
    
    When the high-symmetry points of two slabs are combined by finding all the
    combinations (e.g. ontop_1-hollow_2, ontop_1-bridge_1, ...) by a 
    GetInterfaceHS a lot of duplicates might be created. Here we use a
    pymatgen.analysis.structure_matcher.StructureMatcher to get rid of these
    duplicates both in the unique and the replicated high symmetry points.
    

    Parameters
    ----------
    hs_u : dict
        Unique high-symmetry points of the interface.
    hs_a : dict
        All high-symmetry points of the interface.
    top_slab : pymatgen.core.surface.Slab or pymatgen.core.structure.Structure
        The top slab of the interface
    bot_slab : pymatgen.core.surface.Slab or pymatgen.core.structure.Structure
        The bottom slab of the interfaces
    structure_matcher : pymatgen.analysis.structure_matcher.StructureMatcher
        Class to find equivalent structures (mirrors, rotations, etc...)

    Returns
    -------
    hs_u : dict
        Unique high Symmetry points of the interface without equivalent entries.
    hs_a : dict
        All high Symmetry points of the interface without equivalent entries.

    """
    structure_list = {}
    for key, value in hs_u.items():
        x_shift = value[0][0]
        y_shift = value[0][1]
        inter_struct = StackAlignedSlabs(bot_slab,
                                         top_slab,
                                         top_shift = [x_shift, y_shift, 0])
        clean_struct = CleanUpSiteProperties(inter_struct)
        structure_list[key] = clean_struct
        
    equivalent_structs = {}
    doubles_found = []
    for name, struct in structure_list.items():
        for name_2, struct_2 in structure_list.items():
            if name != name_2:
                if structure_matcher.fit(struct, struct_2) and name not in doubles_found:
                    equivalent_structs.setdefault(name, []).append(name_2)
                    doubles_found.append(name_2)
                    
    for value in equivalent_structs.values():
        for key in value:
            hs_u.pop(key)
            hs_a.pop(key)
            
    return hs_u, hs_a
    

def NormalizeHSDict(hs, to_array=True):
    """
    Convert the hs elements returned by GetSlabHS to proper np.array or lists 
    Important to use with unfolded HS points, which are lists of arrays
    
    """
    
    hs_new = {}
    
    for k in hs.keys():
        data = hs[k]
        
        # Default: convert everything in a list
        if isinstance(data, list):
            elements_list = []
            for element in data:
                elements_list.append(list(element))
            hs_new[k] = elements_list
        else:
            hs_new[k] = [list(data)]
        
        # Convert to array, default and smart choice
        if to_array == True:
            hs_new[k] = np.array(hs_new[k])
    
    return hs_new


def PBC_HSPoints(hs, cell, to_array=False, z_red=True):
    """
    Create a "fake" molecule structure from the HS points calculated for
    the tribological interface, in order to apply PBC to the sites.
    Return a dictionary with the sites within cell. z_red remove the 
    z-coordinates from the translations.

    """
    
    # Type check and error handling
    if not isinstance(hs, dict):
            raise TypeError("hs must be a dictionary")
    
    # Convert the hs elements to lists if necessary
    typ = list( set(type(k) for k in hs.values()) )
    if len(typ) > 1:
        raise ValueError('Your dictionary is weird, values have mixed types')
    elif typ[0] != list:
        hs = HS_DictConverter(hs, to_array=False)
    
    # Run over dictionary values and apply PBC
    hs_new = {}
    for k in hs.keys():
        sites = hs[k]
        
        # Create a fake atomic structures and apply PBC
        atoms_fake = Atoms( positions=sites, cell=cell, pbc=[1,1,1] )
        hs_new[k] = atoms_fake.get_positions( wrap=True, pbc=True )
        
        # Remove z coordinates
        if z_red:
            hs_new[k] = hs_new[k][:, :2]
            
    hs_new = HS_DictConverter(hs_new, to_array=to_array)
    
    return hs_new    
