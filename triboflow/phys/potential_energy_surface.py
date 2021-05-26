#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 18:12:33 2020

Python functions to get the Potential Energy Surface (PES) of an interface.

The module contains the following functions:

    - get_pes
    - remove_duplicates
    - unfold_pes
    - plot_pes

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna
    Code readapted from our past homogeneous workflow, MIT license,
    https://github.com/mcrighi/interface-workflow

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__credits__ = 'Code readapted from our past homogeneous workflow, MIT license, https://github.com/mcrighi/interface-workflow,'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 8th, 2021'

import numpy as np
from scipy.interpolate import Rbf

from triboflow.utils.phys_tools import replicate_points, generate_uniform_grid,\
    orthorombize, pbc_coordinates


# =============================================================================
# EVALUATION OF THE PES - MAIN
# =============================================================================


def get_pes(hs_all, E, cell, to_fig=None, point_density=20):
    """Interpolate the PES using a list of high symmetry points and energies.
    
    Main function to get the Potential Energy Surface (PES) for an interface. 
    The points are replicated to span a 3x3 lattice cell and are interpolated
    by using Radial Basis Functions (cubic function).
    In the output data the energy is normalized so that the absolute minimum
    is 0. Furthermore it is made sure that the lateral point are inside the
    unit cell before they are replicated.
    
    Parameters
    ----------        
    hs : dict
        Unfolded HS points of the interface, covering all the surface.
        
    E : dict
        Contains the energy calculated for each unique surfacial HS site.
    
    cell : numpy.ndarray
        Vectors of the lattice cell of the interface.
        
    #to_fig : string, optional CURRENTLY NOT IMPLEMENTED
    #    Name of the image that you want to save, it will be: 'name_fig'+'.pdf' 
    #    Suggested name: name = 'PES_' + 'Name of the interface'.
    #    The default is None and no image is saved.

    Returns
    -------
    rbf : scipy.interpolate.rbf.Rbf
        Object containing the information of the interpolation of the potential
        energy. Call it on a set of [x, y] coordinates to obtain the energy 
        values for those points. Usage: rbf([x, y])
        
    pes_dict : dict
        Dictionary containing the points and the associated energies divided
        by type. Useful to keep track of the type of the HS sites and their
        energies. To be saved to DB (?)
    
    pes_data : numpy.ndarray
        The entire set of HS points covering the interface with the 
        corresponding energies.
        Format:
            x[0]  y[0]  E[0]
            x[1]  y[1]  E[1]
             .     .     .
             .     .     .
             .     .     .
             
    """
    
    # Unfold the PES points
    E_list, data = unfold_pes(hs_all, E)
    
    #making sure points are not represented twice by ensuring rows in data are unique
    data = remove_duplicates(data)
    #print(len(data))
    
    #make sure that the x and y coordinates are inside the unit cell.
    # x_y_insideCell = pbc_coordinates(data[:, :2],
    #                                  cell,
    #                                  to_array=True)
    # data[:, :2] = x_y_insideCell
        
    # Normalize the minimum to 0
    data[:,2] = (data[:,2]-min(data[:,2]))
    
    # Interpolate the data with Radial Basis Function
    data_rep = replicate_points(data, cell, replicate_of=(3, 3) )
    rbf = Rbf(data_rep[:, 0], data_rep[:, 1], data_rep[:, 2], function='cubic')
    
    # Calculate the PES on a very dense and uniform grid. Useful for further 
    # analysis (MEP, shear strength) and to plot the PES
    coordinates = generate_uniform_grid(cell*2, density=point_density)
    E_new = rbf(coordinates[:, 0], coordinates[:, 1])
    pes_data = np.column_stack([coordinates[:, :2], E_new])
    
    min_x = min(coordinates[:,0])
    max_x = max(coordinates[:,0])
    len_x = max_x - min_x
    dist_x = len_x/(point_density*10)
    min_y = min(coordinates[:,1])
    max_y = max(coordinates[:,1])
    len_y = max_y - min_y
    dist_y = len_y/(point_density*10)
    grid_x = np.arange(min_x, max_x, dist_x)
    grid_y = np.arange(min_y, max_y, dist_y)
    X, Y = np.meshgrid(grid_x, grid_y)
    Z = rbf(X, Y)
    
    to_plot = [X, Y, Z]
    
    return rbf, E_list, pes_data, data, to_plot


# =============================================================================
# UTILITY FOR THE PES
# =============================================================================

def remove_duplicates(data, rounding_decimal=5):
        
    xy = np.round(data[:, :2], decimals=rounding_decimal)
    E_dict = {}
    for i in xy:
        duplicates = np.where((xy == i).all(axis=1))[0]
        E_list = []
        for j in duplicates:
            E_list.append(data[j,2])
        E_dict[str(i)]=np.mean(E_list)
    xy_unique = np.unique(xy, axis=0)
    E=[]
    for i in xy_unique:
        E.append([E_dict[str(i)]])
    
    return np.hstack((xy_unique, np.array(E)))

def unfold_pes(hs_all, E_unique):
    """
    Unfold the energies calculated for the unique HS points of an interface,
    associating them to the replicated HS points covering the whole surface
    cell. hs_all is a dictionary, E_unique a list.

    Parameters
    ----------
        
    hs_all : dict
        Surfacial HS sites that has been unfolded (replicated) across the whole
        lattice cell of slab. 
        Ex. To the key 'ontop_1 + bridge_1' will correspond n points, spanning
        the entire plane axb of the lattice cell. Data is a (n, 3) numpy array.
    
    E : list
        Contains the energy calculated for each unique interfacial HS site.
        The energy are calculated by means of ab initio simulations (VASP).
        E must have the following structure:
            
            [ [label_1, x_1, y_1, E_1], 
              [label_2, x_2, y_2, E_2], 
              ...                      ]
        
        Ex. label should corresponds to the keys in hs_all, associated to a 
        certain shit between the lower and upper slab, e.g. 'ontop_1+bridge_1'.
        
    Returns
    -------
    E_list : list
        It's basically the same as E_unique but contains all the HS points 
        replicated on the whole interface. The structure of the list is:
            
            [ [label_1, x_1, y_1, E_1], 
              [label_2, x_2, y_2, E_2], 
              ...                      ]
        
    E_array : np.ndarray
        Numpy matrix containing the coordinates and the energy useful to 
        interpolate the PES. It's E_list without labels and with array type. 
        The structure of the matrix is:
            
            np.array([ [x_1, y_1, E_1], 
                       [x_2, y_2, E_2], 
                       ...             ])

    """

    # Initialize lists for the result
    E_list = []
    E_array = []
    
    # Extract the element
    for element in E_unique:
       label  = element[0]
       energy = element[3]
       
       # Associate each Energy to all the corresponding HS values
       for row in hs_all[label]:
          x_shift = row[0]
          y_shift = row[1]
          
          E_list.append([label, x_shift, y_shift, energy])
          E_array.append([x_shift, y_shift, energy])
          
    E_array = np.array(E_array)      
    
    return E_list, E_array

def plot_pes(data, lattice, to_fig=None):
    """
    Plot the PES and eventually save it

    """
    
    import matplotlib.pyplot as plt
    
    a = lattice[0]
    b = lattice[1]
    x = data[:, 0]
    y = data[:, 1]
    E = data[:, 2]
    
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


# =============================================================================
# TESTING
# =============================================================================


if __name__ == '__main__':
    print('Testing the creation of a uniform grid for the PES\n')
    vectors = np.array([[3, 0, 0], [0.8, 4, 0.3]])
    a = generate_uniform_grid(vectors, density=1, pts_a=5, to_plot=True)
