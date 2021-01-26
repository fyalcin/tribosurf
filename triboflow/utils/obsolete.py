#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
This file contains obsolete HS, PESfunctions. Might be of some utility to take
inspiration to write other functions or used in rare circumstances

Created on Wed Oct 28 16:20:01 2020

@author: gl
"""

import numpy as np


# =============================================================================
# GENERIC TOOLS - OBSOLETE
# =============================================================================


# Substituted by an conditional statement in PBC_HSPoints
def RemoveZCoords(hs):
    """
    Remove the z coordinates from the HS points saved in a dictionary

    """
    
    for k in hs.keys():
        hs[k] = hs[k][:, :2]
        
    return hs


# It works but it's cumbersome. Understand if a point is inside a cell
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


# Dependency of PointIsInsideLattice
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


# Dependency of CreateBoundary
def IntermediatesPts(p1, p2, npts=100):
    """"
    Return an array of npts equally spaced between p1 and p2.
    Dependency of CreateBoundary
    
    """
    
    delta_x = (p2[0] - p1[0]) / (npts + 1)
    delta_y = (p2[1] - p1[1]) / (npts + 1)
    delta_z = (p2[2] - p1[2]) / (npts + 1)
    
    pts = [ [p1[0] + i*delta_x, p1[1] + i*delta_y, p1[2] + i*delta_z] 
            for i in range(1, npts+1) ]
    
    return np.array(pts)


# Dependency of PointIsInsideLattice
def ClosestPoint(S, q):
    """
    Find the closest point in the set of points S to q.
    Dependency of PointIsInsideLattice
    
    """
    
    closer_pts = np.zeros(3)
    distance = np.sqrt( np.sum((q-S)*(q-S), axis=-1) )
    
    return S[np.argmin(distance), :]


# =============================================================================
# TOOLS TO DEAL WITH PERIODIC BOUNDARY CONDITIONS - OBSOLETE
# =============================================================================


# It works but require nnp which requires pytorch (`700 Mb)
def ApplyPbcToHS(slab, hs): 
    """
    Apply pbc to the HS points of a slab/structure object in the axb plane.
    WARNING: OBSOLETE FUNCTION. It requires nnp and pytorch to work.
    
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


# =============================================================================
# UTILITY FOR THE PES - OBSOLETE
# =============================================================================


def UnfoldPES_FromDict(E, hs_all):
    """
    Unfold the energies calculated for the unique HS points of an interface,
    associating them to the replicated HS points covering the whole surface
    cell. It uses a dictionary for both E and hs_all

    Parameters
    ----------
    E : dict
        Contains the energy calculated for each Hs site of the interface.
        Ex. E may contains the key 'ontop_1 + bridge_1', corresponding to a 
        certain shift between lower and upper slab. The value of the dictionary
        is the ab initio, equilibrium energy obtained for that configuration.
        
    hs_all : dict
        Surfacial HS sites that has been unfolded (replicated) across the whole
        lattice cell of slab. 
        Ex. To the key 'ontop_1 + bridge_1' will correspond n points, spanning
        the entire plane axb of the lattice cell. Data is a (n, 3) numpy array.
        
    Returns
    -------
    pes_data : np.ndarray
        Numpy matrix containing the coordinates and the energy useful to 
        interpolate the PES. The structure of the matrix is:
            
            x[0]  y[0]  E[0]
            x[1]  y[1]  E[1]
             .     .     .
             .     .     .
             .     .     .
        
    pes_dict : dict
        Dictionary containing the coordinates and the corresponding energies
        associated to each HS point of the interface.

    """
    
    pes_dict = {}
    pes_data = []
    
    # WARNING: The elements of hs_all should not have the z coordinates.
    # Call RemoveZCoords() before using this function
    for k in hs_all.keys():
        data = np.column_stack((hs_all[k], np.full(hs_all[k], E[k])))
        pes_dict[k] = data
        pes_data.append(data)
    
    pes_data = np.concatenate(pes_data)
    
    return pes_data, pes_dict