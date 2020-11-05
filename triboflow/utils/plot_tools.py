#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:46:24 2020

Utility tools to calculate the High Simmetry (HS) points for slab and interface

@author: gl
"""

import numpy as np
import matplotlib.pyplot as plt  
from mpl_toolkits.mplot3d import Axes3D
from pymatgen.analysis.adsorption import plot_slab
from triboflow.phys.high_symmetry import HS_DictConverter


# =============================================================================
# PLOTTING TOOLS
# =============================================================================


def Plot_SlabHS(hs, slab, to_fig=None):
    """
    Plot the slab, displaying the atoms and the HS sites of the surface
    
    Parameters
    ----------
    slab : pymatgen.core.surface.Slab 
        The slab to be displayed
        
    hs : dict
        HS sites of the slab.
        
    name_fig : string, optional
        Name of the image that you want to save, it will be: 'name_fig'+'.pdf' 
        Suggested name: 'Name of the material' + 'Miller index'.
        The default is None and no image is saved.     
        
    Returns
    -------
    None.

    """
    
    # Check the type of the hs points
    typ = list( set(type(k) for k in hs.values()) )[0]
    if typ == list:
        hs = HS_DictConverter(hs, to_array=True)
    
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


def Plot_PES(data, lattice, to_fig=None):
    """
    Plot the PES and eventually save it

    """

    a = lattice[0]
    b = lattice[1]
    x = data[0]
    y = data[1]
    E = data[2]

    fact=1.
    level= 43
    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    anglerot='vertical'
    shrin=1.
    zt1=plt.contourf(x, y, E, level, extent=(-fact*a, fact*a, -fact*b, fact*b), cmap=plt.cm.RdYlBu_r)
    cbar1=plt.colorbar(zt1,ax=ax,orientation=anglerot,shrink=shrin)
    cbar1.set_label(r'$E_{adh} (J/m^2)$', rotation=270, labelpad=20,fontsize=15,family='serif')

    ax.quiver(0. , 0., 1., 0.,scale=1.,scale_units='inches',width=0.01,color='white')
    ax.quiver(0. , 0., 0., 1.,scale=1.,scale_units='inches',width=0.01,color='white')
    ax.plot(0.,0.,'w.',ms=7)
    ax.text(0.5,-0.5,'[1 0 1]',rotation='horizontal',color='white', fontsize=14)
    ax.text(-0.5,1.,'[1 2 1]',rotation='vertical',color='white', fontsize=14)
    ax.axis([-fact*a, fact*a, -fact*b, fact*b])
    plt.xlabel(r"distance ($\AA$)",fontsize=12,family='serif')
    plt.ylabel(r"distance ($\AA$)",fontsize=12,family='serif')

    for zt1 in zt1.collections:
       zt1.set_edgecolor("face")
       zt1.set_linewidth(0.000000000001)

    if to_fig != None:
        plt.title("PES for " + str(to_fig), fontsize=18, family='serif')
        plt.savefig('PES_' + str(to_fig) + '.pdf', dpi=300)

    
def Plot_UniformGrid(grid, cell, n_a, n_b):
    """
    Plot an uniform grid of n_aXn_b points on the planar base of a lattice 
    
    """
    
    a = cell[0, :]
    b = cell[1, :]
    v = np.cross(a, b)
    
    mod_a = np.sqrt(a[0]**2. + a[1]**2. + a[2]**2.)
    mod_b = np.sqrt(b[0]**2. + b[1]**2. + b[2]**2.)
    A = np.sqrt(v[0]**2. + v[1]**2. + v[2]**2.)
    
    N = n_a * n_b
    density = N / A
    
    # Print information
    print("1st vector:  {:} -> norm: {:.3f}".format(a, mod_a))
    print("2nd vector:  {:} -> norm: {:.3f}".format(b, mod_b))
    print("N pts: {:}   Area: {:.3f}   Density: {:.3f}".format(N, A, density))
    print("\nUniform {0}x{1} grid\n".format(n_a, n_b))
    print(grid)      
    
    # Projection on the plane, top view
    plt.title("Projection on xy plane")
    plt.plot(grid[:, 0], grid[:, 1], 'o')
    
    # 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(grid[:,0], grid[:,1], grid[:,2], 
               c='r', marker='o', label="3D grid")
    
    # Plot the lattice edge of the plane
    x = [0, a[0], a[0]+b[0], b[0],0]
    y = [0, a[1], a[1]+b[1], b[1],0]
    z = [0, a[2], a[2]+b[2], b[2],0]
    
    ax.plot(x, y, z)
    plt.show()
