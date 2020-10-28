#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Test PES/shear strength modules, using test points
Created on Mon Oct  26 11:20:07 2020

@author: gl
"""

import numpy as np
from PES_functions import ReplicatePESPoints
from utility_functions import PBC_Coordinates


# Define input file containing the data, and the lattice cell vectors
input_file = 'pes_clean.dat'
cell = np.array([[  1.23817 ,   2.1466  ,   0.065971],
                 [  2.478096,   0.      ,   0.065971],
                 [  0.      ,   0.      , -46.575684]])

# Clean the pes data
pes_data = np.genfromtxt('pes_clean.dat')
pes_data = np.unique(pes_data, axis=0)
pbc = PBC_Coordinates(pes_data[:, :2], cell, to_array=True)
pes_data[:, :2] = pbc.copy()

data_rep = ReplicatePESPoints(pes_data, cell, replicate_of=(3, 3))






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
    
    from PES_functions import GetPES, Orthorombize
    from MEP_functions import GetBSMEP, GetMEP
    from SS_functions import GetShearStrength_xy, GetShearStrength
    
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