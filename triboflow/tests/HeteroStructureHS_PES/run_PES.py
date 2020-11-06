#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Test PES/shear strength modules, using test points
Created on Mon Oct  26 11:20:07 2020

@author: gl
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf
from triboflow.utils.plot_tools import Plot_PES
from triboflow.utils.phys_tools import PBC_Coordinates, ReplicatePoints, \
                              GenerateUniformGrid, Orthorombize


# Define input file containing the data, and the lattice cell vectors
input_file = 'pes_clean.dat'
cell = np.array([[  1.23817 ,   2.1466  ,   0.065971],
                 [  2.478096,   0.      ,   0.065971],
                 [  0.      ,   0.      ,  46.575684]])

# Clean the pes data
pes_array = np.genfromtxt('pes_clean.dat')
pes_array = np.unique(pes_array, axis=0)
pbc = PBC_Coordinates(pes_array[:, :2], cell, to_array=True)
pes_array[:, :2] = pbc.copy()

pes_array[:,2] = pes_array[:,2]-min(pes_array[:,2])

data_rep = ReplicatePoints(pes_array, cell, replicate_of=(4, 4))
rbf = Rbf(data_rep[:, 0], data_rep[:, 1], data_rep[:, 2], function='cubic')
coordinates = GenerateUniformGrid(cell, density=30)
x = coordinates[:, 0]
y = coordinates[:, 1]
#E_new = rbf(x, y)
#E_new.reshape(len(x), len(y))

#pes_data = np.column_stack([coordinates[:, :2], E_new])

coords, cell_ortho = Orthorombize(coordinates[:, :2], cell)
#x=coords[:, 0].reshape(49,4)
#y=coords[:, 0].reshape(49,4)
#E=rbf(x, y)#.reshape(len(x), len(y))
#data =  np.array(np.column_stack([x, y, E]))

#Plot_PES(data, cell_ortho, to_fig=None)


orth_coords, cell_ortho = Orthorombize(data_rep, cell)
# ortho_cell_3d = np.zeros((3,3))
# ortho_cell_3d[0,0] = cell_ortho[0,0]
# ortho_cell_3d[1,1] = cell_ortho[1,1]
coordinates = GenerateUniformGrid(cell_ortho, density=30)

#plt.plot(coordinates[:,0], coordinates[:,1], 'ro')
plt.arrow(0,0, cell_ortho[0,0], cell_ortho[0,1])
plt.arrow(0,0, cell_ortho[1,0], cell_ortho[1,1])

X, Y = np.meshgrid(coordinates[:,0], coordinates[:,1])

E = rbf(X, Y)

Plot_PES([X, Y, E], cell_ortho)



# =============================================================================
# MAIN - TO TEST EVERYTHING
# =============================================================================


def StaticTribo(hs, E, cell):  
    """
    Main function to calculate the PES, MEP, and shear strength of an interface

    Parameters
    ----------        
    hs : dict
        Unfolded HS points of the interface, covering all the surface.        

    E : dict
        Contains the energy calculated for each unique surfacial HS site.
        
    cell : numpy.ndarray
        Lattice parameter of the interface you want to study

    Returns
    -------
    pes_info : list
        Relevant information about PES
        
    mep_info : list
        Relevant information about MEP
        
    ss_info : list
        Relevant information about Shear Strength
    """
    
    from triboflow.utils.PES import GetPES, Orthorombize
    from triboflow.utils.MEP import GetBSMEP, GetMEP
    from triboflow.utils.SS import GetShearStrength_xy, GetShearStrength
    
    # Get the PES
    rbf, pes_dict, pes_data = GetPES(hs, E, cell)
    data_ortho, cell_ortho = Orthorombize(pes_data, cell)
    
    # Get the MEP on the potential energy surface starting from a guess
    bsmep, ss_bsmep, theta = GetBSMEP(cell_ortho, rbf)
    mep, mep_convergency = GetMEP(cell_ortho, rbf, theta)
    
    # Calculate the MEP along the x, y, and MEP directions
    data_ss_xy, ss_xy = GetShearStrength_xy(cell_ortho, rbf)
    data_ss_mep, ss_mep  = GetShearStrength(mep, rbf)
    
    pes_info = [pes_dict, pes_data, rbf]
    mep_info = [data_ss_mep, mep_convergency, mep, bsmep]
    ss_info = [ss_xy, ss_bsmep, ss_mep]
    
    return pes_info, mep_info, ss_info


# =============================================================================
# DOCUMENTATION AND COMMENTS
# =============================================================================


# DESCRIPTION OF AN INTERFACE OBJECT, FOR DOCUMENTATION
#interface : pymatgen.core.surface.Slab (pymatgen.core.structure.Structure)
        #The interface object of which you want to calculate the PES. It is 
        #needed to extract the lattice parameter used to unfold the energies.
        #type(slab) could be either a Slab object or a Structure object.