#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 18:12:33 2020

Python functions to get the Potential Energy Surface (PES) of an interface

@author: gl
"""

import numpy as np
from utility_functions import Plot_UniformGrid


# =============================================================================
# EVALUATION OF THE PES - MAIN
# =============================================================================


def GetPES(hs_all, E, cell, to_fig=None):
    """
    Main function to get the Potential Energy Surface (PES) for an interface. 
    The points are replicated to span a 3x3 lattice cell and are interpolated
    by using Radial Basis Functions (cubic function).
    
    ----------        
    hs : dict
        Unfolded HS points of the interface, covering all the surface.
        
    E : dict
        Contains the energy calculated for each unique surfacial HS site.
    
    cell : numpy.ndarray
        Vectors of the lattice cell of the interface.
        
    to_fig : string, optional
        Name of the image that you want to save, it will be: 'name_fig'+'.pdf' 
        Suggested name: name = 'PES_' + 'Name of the interface'.
        The default is None and no image is saved.

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
    
    from scipy.interpolate import Rbf
    
    # Unfold the PES points
    data, pes_dict = UnfoldPES(E, hs_all) 
    
    # Interpolate the data with Radial Basis Function
    data_rep = ReplicatePESPoints(data, cell, replicate_of=(3, 3) )
    rbf = Rbf(data_rep[:, 0], data_rep[:, 1], data_rep[:, 2], function='cubic')
    
    # Calculate the PES on a very dense and uniform grid. Useful for further 
    # analysis (MEP, shear strength) and to plot the PES
    coordinates = GenerateGridForPES(cell, density=10)
    E_new = rbf(coordinates[:, :2])
    pes_data = np.column_stack([coordinates[:, :2], E_new])
    
    return rbf, pes_dict, pes_data


# =============================================================================
# UTILITY FOR THE PES
# =============================================================================


def UnfoldPES(hs_all, E_unique):
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


def ReplicatePESPoints(pes_data, cell, replicate_of=(1, 1)):
    """ 
    Replicate the PES points to cover a (n,m)-size lattice cell
    
    """
    
    n = int(replicate_of[0])
    m = int(replicate_of[1])
    
    # Check wether the number inserted are correct
    if n<=0: n=1
    if m<=0: m=1
    
    if n == 1 and m == 1:
        return pes_data
    
    else:        
        a = cell[0, :]
        b = cell[1, :]
        
        x = pes_data[:, 0]
        y = pes_data[:, 1]
        E = pes_data[:, 2]     
        
        x_new = np.array([])
        y_new = np.array([])
        E_new = np.array([])    
        
        for i in range(n):
                for j in range(m):
                    
                    # Replicate the x- and y- coordinates
                    x_add = x + a[0]*i + b[0]*j
                    y_add = y + a[1]*i + b[1]*j
                    
                    # Collect coordinates and energies
                    x_new = np.append(x_new, x_add)
                    y_new = np.append(y_new, y_add)
                    E_new = np.append(E_new, E)
        
        coordinates_new = np.column_stack([x_new, y_new, E_new])
    
        return coordinates_new


def GenerateGridForPES(cell, density=1, pts_a=None, to_plot=False):
    """
    Generate a 2D-uniform grid of points of density=density on a lattice plane
    given by lattice[0,:]Xlattice[1,:]

    Parameters
    ----------
    lattice : numpy.ndarray
        Vectors of the lattice cell. A uniform grid of points is generated on
        the surface spanned by the first and second vector, i.e. axb.
        lattice shape is (2, 3) or (3, 3), the third vector is not necessary.
        lattice is in Angstrom units.

    density : int, optional
        Density of the grid of points that will cover the planar surface of 
        the lattice cell. Units: number of points per unit Angstrom^2
        
    pts_a : int, optional
        If this value is provided, the grid will contain pts_a points along 
        the first vector and (b/a)*pts_a along the second vector. 
        a,b : lengths of the planar lattice vectors. The default is None.
                
    to_plot : bool, optional
        Wether to display the grid of points inside the lattice cell. 
        Plot is redirected to standard output. The default is False.

    Returns
    -------
    matrix : numpy.ndarray
        Grid of points spanning the entire lattice plane.
        Format:
            x[0]  y[0]  z[0]
            y[1]  y[1]  z[1]
             .     .     .
             .     .     .
             .     .     .

    """
        
    a = cell[0, :]
    b = cell[1, :]
    a_mod = np.sqrt(a[0]**2. + a[1]**2. + a[2]**2.)
    b_mod = np.sqrt(b[0]**2. + b[1]**2. + b[2]**2.)
    ratio = b_mod/a_mod
    
    # Calculate the number of points for each lattice vector
    if pts_a == None:
        N_tot = round(density * a_mod * b_mod)
        n_a = int(round( np.sqrt( N_tot/ratio )))
        n_b = int(round( ratio*n_a ))
    else:
        n_a = pts_a
        n_b = int(round( ratio*n_a ))
    
    # Obtain the displacements along a and b
    dist_a_x = a[0]/n_a 
    dist_a_y = a[1]/n_a
    dist_a_z = a[2]/n_a
    dist_b_x = b[0]/n_b
    dist_b_y = b[1]/n_b
    dist_b_z = b[2]/n_b
    
    # Create the grid
    matrix = np.zeros((n_a*n_b, 3))
    k = 0
    for i in range(0, n_a):
        for j in range(0, n_b):
            matrix[k, 0] = i*dist_a_x + j*dist_b_x
            matrix[k, 1] = i*dist_a_y + j*dist_b_y
            matrix[k, 2] = i*dist_a_z + j*dist_b_z
            k += 1
    if to_plot:
        Plot_UniformGrid(matrix, cell, n_a, n_b)

    return matrix


def Orthorombize(pes_data, cell):
    """
    Take the replicated points of the pes and cut them in a squared shape.
    TODO : Improve the code and VECTORIZE

    """
    
    a = cell[0, :]
    b = cell[1, :]
    
    if np.sign(a[0]) == np.sign(b[0]):
        if a[0] > 0:
            x_up = a[0] + b [0]
            x_dw = 0
        else:
            x_up = 0
            x_dw = a[0] + b [0]           
    else:
        x_up =  max(abs(a[0]), abs(b[0]))
        x_dw =  min(abs(a[0]), abs(b[0]))
    
    if np.sign(a[1]) == np.sign(b[1]):
        if a[1] > 0:
            y_up = a[1] + b[1]
            y_dw = 0
        else:
            y_up = 0
            y_dw = a[1] 
    else:
        y_up =  max(abs(a[1]), abs(b[1]))
        y_dw =  min(abs(a[1]), abs(b[1]))
        
    index_x = pes_data[:, 0] <= 2*x_up and pes_data[:, 0] >= 2*x_dw
    index_y = pes_data[:, 1] <= 2*y_up and pes_data[:, 1] >= 2*y_dw 
    index = index_x * index_y
    
    orthorombic =  []
    for i, row in enumerate(pes_data):
        if index[i] == True:
            orthorombic.append(row)
    
    orthorombic = np.column_stack(orthorombic)
    cell = np.array([[x_up, y_dw], [x_dw, y_up]])
    
    return orthorombic, cell


# =============================================================================
# TESTING
# =============================================================================


if __name__ == '__main__':
    print('Testing the creation of a uniform grid for the PES\n')
    vectors = np.array([[3, 0, 0], [0.8, 4, 0.3]])
    a= GenerateGridForPES(vectors, density=1, pts_a=5, to_plot=True)