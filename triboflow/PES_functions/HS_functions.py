#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:52:15 2020

Calculate the High Simmetry (HS) points for slab and interface

@author: gl
"""

import numpy as np


# =============================================================================
# CALCULATE THE HS POINTS FOR A SLAB
# =============================================================================

def GetSlabHS(slab, allowed_sites=['ontop', 'bridge', 'hollow']): 
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
        
    normalize : bool, optional
        The default is True.

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

    from pymatgen.analysis.adsorption import AdsorbateSiteFinder
    
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
    hs = NormalizeHSDict(hs, to_array=True)
    hs_all = NormalizeHSDict(hs_all, to_array=True)

    return hs, hs_all


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


# =============================================================================
# CALCULATE HS POINTS FOR AN INTERFACE
# =============================================================================

def GetInterfaceHS(hs_1, hs_2):
    """/fs/home/wolloch/git_test/TriboFlow/triboflow/
    Calculate the HS sites for a hetero interface by combining the HS sites of
    the bottom slab (hs_1) with the upper slab (hs_2) 

    Parameters
    ----------
    hs_1 : dict
        High Symmetry points of the bottom slab
    hs_2 : dict
        High Symmetry points of the upper slab

    Returns
    -------
    hs : dict
        High Symmetry points of the hetero interface

    """
    
    hs = {}
    
    # Calculate the shift between each HS point of the first material with each
    # HS point of the second material
    for k1 in hs_1.keys():
        for k2 in hs_2.keys():          
            shifts_stack = []
            d1 = hs_1[k1]
            d2 = hs_2[k2]
            
            for el_d1 in d1:
                shifts_stack.append( d2 - el_d1 )
                
            hs[k1+'-'+k2] = np.concatenate(shifts_stack, axis=0)
            
    return hs