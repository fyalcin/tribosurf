#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:46:24 2020

Utility tools to calculate the High Simmetry (HS) points for slab and interface

@author: gl
"""

import numpy as np
from HS_functions import *
from ase import Atoms


# =============================================================================
# GENERIC TOOLS
# =============================================================================


def RemoveZCoords(hs):
    """
    Remove the z coordinates from the HS points saved in a dictionary

    """
    
    for k in hs.keys():
        hs[k] = hs[k][:, :2]
        
    return hs


# =============================================================================
# PLOT TOOLS
# =============================================================================


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
        The default is None and no image is saved.
         
        
    Returns
    -------
    None.

    """
 
    import matplotlib.pyplot as plt  
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
    
    
def Plot_UniformGrid(lattice, grid, n_a, n_b):
    """
    Plot an uniform grid of n_aXn_b points on the planar base of a lattice 
    
    """
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    a = lattice[0, :]
    b = lattice[1, :]
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


# =============================================================================
# TOOLS TO DEAL WITH PERIODIC BOUNDARY CONDITIONS
# =============================================================================


def ApplyPbcToHS(slab, hs):
    """
    Apply pbc to the HS points of a slab/structure object in the axb plane
    
    """
    
    from torch import Tensor
    import nnp.pbc as pbc
    
    lattice_t = Tensor(slab.lattice.matrix)
    pbc_t     = Tensor([True, True, False])
    hs_new = {}
    
    for k in hs.keys():
        element_new = pbc.map2central( lattice_t, 
                                       Tensor(hs[k]),
                                       pbc_t )
        hs_new[k] = np.array(element_new)
            
    return hs_new

    
def PointIsInsideLattice(lattice, q):
    
    """
    Determine wether a point q is inside a lattice unit cell.
    Return the reponse (True or False), the closest point to q along the closed
    contour made by axb, and the vector at which q is closest (a or b).
    
    TODO:
        - For now only check for the planar quadrilateral given by axb.
          To be generalized to a 3D lattice cell.
        - VECTORIZE the function, to check simultaneously for an array q.
    
    """
    
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
        alpha = ( (q[0]-belong_p[0][0])*(belong_p[1][0]-belong_p[0][0])   +   \
                  (q[1]-belong_p[0][1])*(belong_p[1][1]-belong_p[0][1])   +   \
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
    """
    Given a lattice cell return an array containing all the points along the
    closed contours formed by vector a and b (axb quadrilateral)
    
    """
    
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
    Return an array of npts equally spaced between p1 and p2.
    Dependenceof CreateBoundary
    
    """
    
    delta_x = (p2[0] - p1[0]) / (npts + 1)
    delta_y = (p2[1] - p1[1]) / (npts + 1)
    delta_z = (p2[2] - p1[2]) / (npts + 1)
    
    pts = [ [p1[0] + i*delta_x, p1[1] + i*delta_y, p1[2] + i*delta_z] 
            for i in range(1, npts+1) ]
    
    return np.array(pts)


def ClosestPoint(S, q):
    """
    Find the closest point in the set of points S to q.
    Dependence of PointIsInsideLattice
    
    """
    
    closer_pts = np.zeros(3)
    distance = np.sqrt( np.sum((q-S)*(q-S), axis=-1) )
    
    return S[np.argmin(distance), :]

def ToList(hs, key):
    """
    Take a dictionary as input and key. Convert the element of the 
    dictionary from a (n, 3) array into a list of lists.

    """
    
    my_list = []
    for row in hs[key]:
        my_list.append(list(row))
    return my_list
    
    
def PBCPoints(hs, cell, z_red=True):
    """
    Create a "fake" molecule structure from the HS points calculated for
    the tribological interface, in order to apply PBC to the sites.
    Return a dictionary with the sites within cell. z_red remove the 
    z-coordinates from the translations.

    """
    
    hs_new = {}
    for k in hs.keys():
        sites = ToList(hs, k)
        atoms_fake = Atoms(positions=sites, cell=cell, pbc=[1,1,1])
        hs_new[k] = atoms_fake.get_positions( wrap=True, pbc=True)
    
    if z_red:
        for k in hs_new.keys():
            hs_new[k] = hs_new[k][:, :2]
    
    return hs_new