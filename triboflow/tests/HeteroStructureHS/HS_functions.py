#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:52:15 2020

Calculate the High Simmetry (HS) points for slab and interface + utility

@author: gl
"""

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
    
    # Add a key to the dictionary related to all the sites
    #hs['all'] = unique_sites['all']
    #hs_all['all'] = replica_sites['all']
    
    # Convert the elements of the dictionaries to proper numpy arrays and 
    # remove the z coordinate of the HS sites
    hs = NormalizeHSDict(hs, to_array=True)
    hs_all = NormalizeHSDict(hs_all, to_array=True)

    return hs, hs_all


# =============================================================================
# CALCULATE HS POINTS FOR AN INTERFACE
# =============================================================================

def GetInterfaceHS(hs_1, hs_2):
    """
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
    # Loop over the keys
    for k1 in hs_1.keys():
        for k2 in hs_2.keys():
            shifts_stack = []
            
            # Check the shape of entry arrays
            d1 = hs_1[k1].reshape((1,3)) if len(hs_2[k2].shape) == 1 \
                                         else hs_1[k1]
            d2 = hs_2[k2].reshape((1,3)) if len(hs_2[k2].shape) == 1 \
                                         else hs_2[k2]
            
            for el_d1 in d1:
                shift = (d2 - el_d1)[:, :2]
                if shift.shape[0] == 1:
                    shift = shift.reshape((2,))
                shifts_stack.append(shift)
            
            hs[k1+' - '+k2] = np.concatenate(shifts_stack, axis=0)
            
    return hs


# =============================================================================
# EVALUATION OF THE PES ON HS POINTS AND UNFOLDING
# =============================================================================


def UnfoldHSPES(E, hs_all):
    """
    Unfold the energies calculated for the unique HS points of an interface,
    by associating them to the replicated HS points of the interface of the
    lattice cell.

    Parameters
    ----------
    E : dict
        Contains the energy calculated for each surfacial HS site.
        
    hs_all : dict
        Contain the surfacial HS sites of the interface replicated to the whole
        lattice cell.

    Returns
    -------
    pes_dict : dict
        Energy grid to compute the PES

    """
    
    pes_all = {}
    
    # Find the shifts to be done
    for k in hs_all.keys():
        pes_all[k] = np.column_stack((hs_all[k], np.full(hs_all[k], E[k])))    
    
    return pes_all


# =============================================================================
# UTILITIES AND PLOT TOOLS
# =============================================================================


def NormalizeHSDict(hs, to_array=True):
    """
    Convert the hs elements returned by GetSlabHS to proper np.array or lists 
    Relevant for unfolded HS sites, which are returned as lists of np.arrays.
    
    """
    
    hs_new = {}
    
    for k in hs.keys():
        data = hs[k]
        
        # Default: convert everything in a list
        if isinstance(data, list):
            elements_list = []
            for element in data:
                elements_list.append(list(element[:2]))
            hs_new[k] = elements_list
        else:
            hs_new[k] = list(data[:2])
        
        # Convert to array
        if to_array == True:
            hs_new[k] = np.array(hs_new[k])
    
    return hs_new
    

def Plot_SlabHS(slab, hs, to_fig=None):
    """
    Plot the slab slab, displaying the atoms and the HS sites of the surface
    
    Parameters
    ----------
    slab : pymatgen.core.surface.Slab 
        The slab to be displayed
        
    hs : dict
        HS sites of the slab.
        
    name_fig : string, optional
        Name of the image that you want to save, it will be: 'name_fig'+'.pdf' 
        Suggested name: 'Name of the material' + 'Miller index'.
        The default is None and no images are saved.
         
        
    Returns
    -------
    None.

    """
 
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from mpl_toolkits.mplot3d import Axes3D    
    from pymatgen.analysis.adsorption import plot_slab
    
    # Extract the lattice vector of the basal plane
    a = slab.lattice.matrix[0, :]
    b = slab.lattice.matrix[1, :]
    
    # plot the atoms and the lattice cell
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plot_slab(slab, ax, scale=0.8, repeat=3, window=1.25, 
              draw_unit_cell=True, decay=0.2, adsorption_sites=False)
    ax.set(xlim = ( -0.1*(a[0] + b[0]), 1.1*(a[0] + b[0]) ), 
           ylim = ( -0.1*(a[1] + b[1]), 1.1*(a[1] + b[1]) ))
    
    # Add the HS sites with the proper labels
    for k in hs.keys():
        data = hs[k]
        if len(data.shape) == 1:
            plt.plot(data[0], data[1], marker='x', markersize=10, mew=3, 
                     linestyle='', zorder=10000, label=k)     
        else:
            plt.plot(data[:,0], data[:,1], marker='x', markersize=7.5, mew=3, 
                     linestyle='', zorder=10000, label=k)
 
    plt.legend(bbox_to_anchor=(1.025, 1), loc='upper left')
    
    if to_fig != None:
        plt.savefig(to_fig+'.pdf', dpi=300)
    
    plt.show()
    
    
def PointIsInsideLattice(lattice, q):
    
    a = lattice[0, :]
    b = lattice[1, :]
    q = np.array(q)
    
    # Calculate the lattice boundary and find the closest point to q
    S, lines = CreateBoundary(lattice)
    p = ClosestPoint(S, q)   
    
    r = q - p
    
    if np.sqrt(r.dot(r)) < 1e-10:
        return True # The point is inside!
    else:
        for n, line in enumerate(lines):
            if (line - p == [0, 0, 0]).all() == 0:
                if n == 0:
                    belong_p = [np.zeros(3), a]
                    vector = a
                elif n == 1:
                    belong_p = [a, a+b]
                    vector = b
                elif n == 2:
                    belong_p = [a+b, b]
                    vector = a
                else:
                    belong_p = [b, np.zeros(3)]
                    vector = b
                break
       
        # Hypothesis to find out the normal to the vector where p belongs, 
        # and passing through q
        alpha = ( (q[0]-belong_p[0][0])*(belong_p[1][0]-belong_p[0][0])   +  \
                  (q[1]-belong_p[0][1])*(belong_p[1][1]-belong_p[0][1])   +  \
                  (q[2]-belong_p[0][2])*(belong_p[1][2]-belong_p[0][2]) )     \
                   /                                                          \
        ( (belong_p[1][0]-belong_p[0][0])*(belong_p[1][0]-belong_p[0][0]) +   \
          (belong_p[1][1]-belong_p[0][1])*(belong_p[1][1]-belong_p[0][1]) +   \
          (belong_p[1][2]-belong_p[0][2])*(belong_p[1][2]-belong_p[0][2]) )

        intersect = [ belong_p[0][0]+alpha*(belong_p[1][0]-belong_p[0][0]), 
                      belong_p[0][1]+alpha*(belong_p[1][1]-belong_p[0][1]),
                      belong_p[0][2]+alpha*(belong_p[1][2]-belong_p[0][2]) ]
        intersect = np.array(intersect)
        
        normal = np.array(intersect) / np.sqrt(intersect.dot(intersect))
    
        # Check the scalar product between the normal to the surface and r
        if np.dot(r, normal) >= 0:
            isinside=False
        else:
            isinside=True
    return not (np.dot(r, normal) >= 0), p, vector

def CreateBoundary(lattice, step=0.05):
    
    a = lattice[0, :]
    b = lattice[1, :]
    n_a = int (np.sqrt(a.dot(a)) / step)
    n_b = int (np.sqrt(b.dot(b)) / step)
    
    # Create the boundaries
    line1 = IntermediatesPts(np.zeros(3), a, n_a)
    line2 = IntermediatesPts(a, a+b, n_b)
    line3 = IntermediatesPts(a+b, b, n_a)
    line4 = IntermediatesPts(b, np.zeros(3), n_b)
    
    S = np.concatenate((np.zeros((1,3)), line1, [a],
                        line2, [a+b], 
                        line3, [b],
                        line4))
    
    return S, [line1, line2, line3, line4]


def IntermediatesPts(p1, p2, npts=100):
    """"
    Return an array of npts equally spaced between p1 and p2
    
    """
    
    delta_x = (p2[0] - p1[0]) / (npts + 1)
    delta_y = (p2[1] - p1[1]) / (npts + 1)
    delta_z = (p2[2] - p1[2]) / (npts + 1)
    
    pts = [ [p1[0] + i*delta_x, p1[1] + i*delta_y, p1[2] + i*delta_z] 
            for i in range(1, npts+1) ]
    
    return np.array(pts)


def ClosestPoint(S, q):
    """
    Find the closest point in the set of points S to q
    
    """
    
    closer_pts = np.zeros(3)
    distance = np.sqrt( np.sum((q-S)*(q-S), axis=-1) )
    
    return S[np.argmin(distance), :]
    
    