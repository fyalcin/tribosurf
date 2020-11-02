#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:46:24 2020

Utility tools to calculate the High Simmetry (HS) points for slab and interface

@author: gl
"""

import numpy as np
from triboflow.PES_functions.HS_functions import HS_DictConverter


# =============================================================================
# TOOLS FOR PERIODIC BOUNDARY CONDITIONS
# =============================================================================

    
def PBC_Coordinates(data, cell, to_array=True, scaled_positions=False):
    """
    Apply Periodic Boundary Conditions to a set of atoms (data) in a given 
    lattice cell (cell). PBC are applied by creating a "fake" molecule.
    Return a numpy array containing the atomic sites within cell.

    """
    
    from ase import Atoms
    
    # Check the types of the input parameters
    if not ( (isinstance(data, list)) or (isinstance(data, np.ndarray)) ):
            raise TypeError("data must be a numpy array or a list")
    
    if not ( (isinstance(cell, list)) or (isinstance(cell, np.ndarray)) ):
            raise TypeError("cell must be a numpy array or a list")
    
    # Manage the situation where you provide only x, y coordinates
    data = np.array(data)
    two_col = False
    pbc = [1, 1, 1]
    if data.shape[1] == 2:
        two_col = True
        pbc = [1, 1, 0]
        data = np.column_stack( [data, np.zeros(len(data))] )
    
    # Create a fake atomic structures and apply PBC
    if scaled_positions:
        atoms_fake = Atoms( scaled_positions=data, cell=cell, pbc=pbc )
    else:
        atoms_fake = Atoms( positions=data, cell=cell, pbc=pbc )
        
    # Get the coordinates
    data_new = atoms_fake.get_positions( wrap=True, pbc=True )
    if two_col:
        data_new = data_new[:, :2]
    
    if not to_array:
        data_new = data_new.tolist()  

    return data_new  


# =============================================================================
# TOOLS TO REPLICATE POINTS IN THE CELL
# =============================================================================


def ReplicatePoints(data, cell, replicate_of=(1, 1)):
    """ 
    Replicate a set of points or atomic sites in a (n,m)-size lattice cell.
    Useful to replicate data to interpolate the PES.
    
    """
    
    n = int(replicate_of[0])
    m = int(replicate_of[1])
    
    # Check wether the number inserted are correct
    if n<=0: n=1
    if m<=0: m=1
    
    if n == 1 and m == 1:
        return data
    
    else:        
        a = cell[0, :]
        b = cell[1, :]
        
        x = data[:, 0]
        y = data[:, 1]
        E = data[:, 2]     
        
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


def GenerateUniformGrid(cell, density=1, pts_a=None, to_plot=False):
    """
    Generate a 2D-uniform grid of points of points with a given density on the
    basal lattice plane of a cell (a X b), i.e. lattice[0,:] X lattice[1,:].
    You can set a uniform density or provide the points along a.

    Parameters
    ----------
    lattice : numpy.ndarray
        Vectors of the lattice cell. A uniform grid of points is generated on
        the surface spanned by the first and second vector, i.e. a X b.
        lattice shape is (2, 3) or (3, 3); lattice is in Angstrom units.

    density : float, optional
        Density of the grid of points that will cover the planar surface of 
        the lattice cell. Units: number of points per unit Angstrom^2
        
    pts_a : int, optional
        If this value is provided, the grid will contain pts_a points along 
        the first vector and (b/a)*pts_a along the second vector. a and b are 
        the lengths of the planar lattice vectors. The default is None.
                
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


# =============================================================================
# TOOLS TO MODIFY THE LATTICE CELL
# =============================================================================


def Orthorombize(data, cell):
    """
    Take the replicated points of the pes and cut them in a squared shape.
    TODO : Improve the code and VECTORIZE

    """
    
    a = cell[0, :]
    b = cell[1, :]
    
    if np.sign(a[0]) == np.sign(b[0]):
        if a[0] > 0:
            x_up = a[0] + b[0]
            x_dw = 0
        else:
            x_up = 0
            x_dw = a[0] + b[0]           
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
        
    index_x = (data[:, 0] <= 2*x_up) * (data[:, 0] >= 2*x_dw)
    index_y = (data[:, 1] <= 2*y_up) * (data[:, 1] >= 2*y_dw) 
    index = index_x * index_y
    
    orthorombic =  []
    for i, row in enumerate(data):
        if index[i] == True:
            orthorombic.append(list(row))
            
    #orthorombic = np.column_stack(orthorombic)
    orthorombic = np.array(orthorombic)
    cell = np.array([[x_up, y_dw], [x_dw, y_up]])
    
    return orthorombic, cell


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
 
    import matplotlib.pyplot as plt  
    from pymatgen.analysis.adsorption import plot_slab
    
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
    
    
def Plot_UniformGrid(grid, cell, n_a, n_b):
    """
    Plot an uniform grid of n_aXn_b points on the planar base of a lattice 
    
    """
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
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
