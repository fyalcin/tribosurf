#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:52:15 2020

Calculate the High Simmetry (HS) points for slab and interface + utility

@author: gl
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, Slab
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.transformations.standard_transformations import RotationTransformation

from mpinterfaces.transformations import get_aligned_lattices, get_interface


# =============================================================================
# CALCULATE THE HS POINTS FOR A SLAB
# =============================================================================

def GetSlabHS(material, near_reduce=0.01, allowed_sites=['ontop', 'bridge', 'hollow']): 
    """
    Extract the High Simmetry (HS) points for a given slab. It returns the 
    unique HS points for the slab and the unfolded points inside the lattice
    cell

    Parameters
    ----------
    material : pymatgen.core.surface.Slab
        The slab of which you want to calculate the HS points of the surface
    near_reduce : float, optional
        Near reduction threshold to be feed into pymatgen. Default is 0.01
        
    Returns
    -------
    hs : dict
        The unique HS points for the selected slab
    hs_all : dict
        HS points for the selected slab replicated inside the cell
        
    """
    
    adsf = AdsorbateSiteFinder(material)
    
    # Extract the unique HS points for the given surface
    unique_sites = adsf.find_adsorption_sites( distance = 0, 
                                               symm_reduce = 0.01, 
                                               near_reduce = near_reduce, 
                                               no_obtuse_hollow = True )
    
    #Extract all the HS points of the lattice cell for the given surface
    replica_sites = adsf.find_adsorption_sites( distance = 0, 
                                                symm_reduce = 0, 
                                                near_reduce = near_reduce, 
                                                no_obtuse_hollow = True )
    
    # Identify the unique HS points of the slab and unfold them in the cell
    hs = {} 
    for key in allowed_sites:
        if unique_sites[key] != []:
            n=1
            for data in unique_sites[key]:
                hs[key+'_'+str(n)] = data
                n += 1
            
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
    
    #hs['all'] = unique_sites['all']
    #hs_all['all'] = replica_sites['all']
    
    hs = NormalizeHSToList(hs)

    return hs, hs_all


def NormalizeHSToList(hs):
    
    keys = list(hs.keys())
    hs_new = hs.copy()
    
    for k in keys:
        data = hs[k]
        lista = []
        for i in np.matrix(data):
            lista.append(list(data))
        hs_new[k] = np.matrix(lista)


# =============================================================================
# CALCULATE HS POINTS FOR AN INTERFACE
# =============================================================================

def GetInterfaceHS(hs_bot, hs_top):
    
    # Extract the keys from the HS dictionaries
    kbot = list(hs_bot.keys())
    ktop = list(hs_top.keys())
    if 'all' in kbot:
        kbot.remove('all')
    if 'all' in ktop:
        ktop.remove('all')
    
    hs = {}
    for k1 in kbot:
        for k2 in ktop:
            shift_list = []
            
            for d1 in np.matrix(list(hs_bot[k1])):
                for d2 in np.matrix(list(hs_top[k2])):                  
                    shift = np.array(d1 - d2)                      
                    shift_list.append([shift[0,0], shift[0,1]])
             
            hs[k1+' + '+k2] = np.array(shift_list)
            
    return hs



def UnfoldHS(hs, hs_all):
    """
    Unfold the PES Energies, calculated for the unique HS points of the 
    interface (reduced by simmetry), to the whole set of HS points of the 
    super lattice cell adopted.

    Parameters
    ----------
    E : TYPE
        DESCRIPTION.
    hs : TYPE
        DESCRIPTION.
    hs_all : TYPE
        DESCRIPTION.

    Returns
    -------
    pes_dict : TYPE
        DESCRIPTION.

    """
    
    coordinates = []
    
    for k in hs_all.keys():
        unique = hs[k]
        coordinates.append(unique)
        
        for shift in hs_all[k]:
            
            c = shift - unique
            if c != [0, 0]:
                coordinates.append(shift - unique)       
    
    return coordinates


def unfold PES(E, hs):
    
    pes_unfolded = {}
    for k in E.keys():
        


# =============================================================================
# PLOT TOOLS AND UTILITY
# =============================================================================


def Plot_SlabHS(material, hs, to_save=False, name_fig='material'):
    """
    Plot a slab, displaying the HS points on top of the surface
    
    Parameters
    ----------
    material : pymatgen.core.surface.Slab
        The slab of the material of interest
    hs : dict
        The high simmetry point for the material
    to_save : bool, optional
        Wether to save an image of the plot. The default is False
    name_fig : string, optional
        Name of the image to save. The default is 'material'.
        Suggested: 'Name of the material' + 'Miller index'

    Returns
    -------
    None.

    """
 
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from mpl_toolkits.mplot3d import Axes3D    
    from pymatgen.analysis.adsorption import plot_slab
    
    # Extract the HS points from the hs dictionary
    hs = np.array(hs['all'])
    
    # Extract the lattice vector of the basal plane
    a = material.lattice.matrix[0, :]
    b = material.lattice.matrix[1, :]
    
    # Make the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_slab(material, ax, repeat=5, draw_unit_cell=True, adsorption_sites=True)
    plt.plot(hs[:,0], hs[:,1], 'o')
    ax.set(xlim = (-np.sign(a[0])*0.5, a[0]+np.sign(a[0])*0.5), 
           ylim = (-np.sign(b[1])*0.5, b[1]+np.sign(b[1])*0.5))
    
    if to_save:
        plt.savefig(name_fig+'.pdf', dpi=300)
    
        
if __name__ == "__main__":
    pass
    
