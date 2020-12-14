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

def CleanUpHSDicts(hs_unique, hs_all, top_aligned, bottom_aligned,
                   decimals=4):
    c_u = hs_unique.copy()
    c_a = hs_all.copy()
    c_hsp_u, c_hsp_a = RemoveDuplicatesFromHSDicts(c_u,
                                                   c_a,
                                                   decimals=decimals)
    c_hsp_u_reduced, c_hsp_a_reduced = RemoveEquivalentShifts(c_hsp_u,
                                                          c_hsp_a,
                                                          top_aligned,
                                                          bottom_aligned)
    hs_a_out = AssigneAllPoints(c_hsp_u_reduced, 
                                c_hsp_a_reduced,
                                top_aligned,
                                bottom_aligned)
    hs_u_out = c_hsp_u_reduced
    
    return hs_u_out, hs_a_out

def AssigneAllPoints(hs_unique, hs_all, top_aligned, bottom_aligned):

    all_shifts = []
    for key, value in hs_all.items():
        if all_shifts == []:
            all_shifts = value
        else:
            all_shifts = np.concatenate([all_shifts, value], axis=0).tolist()
    all_shifts = np.unique(all_shifts, axis=0)

    struct_match = StructureMatcher(ltol=0.01, stol=0.01, angle_tol=0.01,
                                primitive_cell=False, scale=False)
    top_slab, bot_slab = ReCenterAlignedSlabs(top_aligned,
                                              bottom_aligned,
                                              d=4.5)
    new_hsp_dict_a = {}
    for key, value in hs_unique.items():
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
            if struct_match.fit(unique_struct, test_struct):
                new_hsp_dict_a.setdefault(key, []).append(shift)
    
    return new_hsp_dict_a

def RemoveEquivalentShifts(hs_unique, hs_all, top_slab, bot_slab,
                           ltol=0.01, stol=0.01, angle_tol=0.01,
                           primitive_cell=False, scale=False):
    hs_u = hs_unique.copy()
    hs_a = hs_all.copy()
    top_slab, bot_slab = ReCenterAlignedSlabs(top_slab, bot_slab, d=4.5)
    struct_match = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol,
                                primitive_cell=primitive_cell, scale=scale)
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
                if struct_match.fit(struct, struct_2) and name not in doubles_found:
                    equivalent_structs.setdefault(name, []).append(name_2)
                    doubles_found.append(name_2)
                    
    for value in equivalent_structs.values():
        for key in value:
            hs_u.pop(key)
            hs_a.pop(key)
            
    
    
    return hs_u, hs_a

def RoundPosInDict(high_symm_dict, decimals=5):
    """Round coordinates for high-symmetry points in dict to selected precision.
    
    Convert the coordinate lists in a high-symmetry point dictionary as
    produced by GetSlabHS and GetInterfaceHS to np.arrays and use np.round
    to round the coordinates. Convert back to lists and return.

    Parameters
    ----------
    high_symm_dict : dict
        Dictionary of high-symmetry points for an interface. As computed
        by GetInterfaceHS or GetSlabHS.
    decimals : int, optional
        Selects the number of decimal points the coordinates are rounded to.
        The default is 5.

    Returns
    -------
    dict
        Just as the input dictionary, but with the coordinates rounded to the
        desired precision.

    """
    hs = high_symm_dict.copy()
    new_hs = {}
    for key, value in hs.items():
        pos = np.array(value)
        rounded_array = np.round(np.array(pos), decimals=decimals)
        new_hs[key] = rounded_array
    
    return jsanitize(new_hs)

def RemoveDuplicatesFromHSDicts(hs_unique, hs_all, decimals=5):
    """Remove the duplicates from dictionaries of high-symmetry points.
    
    Sometimes duplicate high-symmetry points with distinct labels are reported
    by GetSlabHS and GetInterfaceHS. This function removes them, while
    ensuring that the identification of high-symmetry points with their label
    is kept correctly and no correct duplicate (points which have the same
    symmetry but different xy-coordinates) are thrown away. The procedure
    also involves rounding the coordinates to the chosen accuracy, so also
    points that are not 'exactly' identical can be identified as equivalent.
    

    Parameters
    ----------
    hs_unique : dict
        Dictionary of unique high-symmetry points for an interface. As computed
        by GetInterfaceHS.
    hs_all : dict
        Dictionary of correctly replicated high-symmetry points for an
        interface. As computed by GetInterfaceHS.
    decimals : int, optional
        Selects the number of decimal points the coordinates are rounded to.
        The default is 5.

    Returns
    -------
    u_hs : dict
        Dictionary of unique high-symmetry points for an interface. Just as the
        input dictionary 'hs_unique', but with wrong duplicates removed.
    a_hs : TYPE
        Dictionary of all high-symmetry points for an interface. Just as the
        input dictionary 'hs_all', but with wrong duplicates removed.

    """
    u_hs = RoundPosInDict(hs_unique.copy(), decimals=decimals)
    a_hs = RoundPosInDict(hs_all.copy(), decimals=decimals)
    rev_dict = {}
    for key, value in u_hs.items():
        #have to make sure that -0.0 is changed to 0.0 to ensure that 
        #equivalent points are recognised...
        value_clean = []
        for i in value[0]:
            if i == 0.0:
                value_clean.append(abs(i))
            else:
                value_clean.append(i)               
        rev_dict.setdefault(str([value_clean]), list()).append(key)
    for value in rev_dict.values():
        if len(value) > 1:
            equivalent_points = np.array(a_hs[value[0]])
            for i in value[1:]:
                #remove the equivalent entries from the u_hs dictionary
                u_hs.pop(i)
                #equivalent_points = np.vstack((equivalent_points,
                #                              np.array(a_hs[i])))
                a_hs.pop(i)
            #equivalent_points = np.unique(equivalent_points, axis=0)
            #a_hs[value[0]] = jsanitize(equivalent_points)
            
    return u_hs, a_hs
    

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


