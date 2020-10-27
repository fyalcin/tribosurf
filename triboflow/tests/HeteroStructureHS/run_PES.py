#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Test PES/shear strength modules, using test points
Created on Mon Oct  26 11:20:07 2020

@author: gl
"""

import numpy as np
from PES_functions import *
from utility_functions import *

pes_array2 = np.genfromtxt('pes_clean.dat')
pes_array = np.unique(pes_array2, axis=0)
cell = np.array([[  1.23817 ,   2.1466  ,   0.065971],
                 [  2.478096,   0.      ,   0.065971],
                 [  0.      ,   0.      , -46.575684]])


data_rep = ReplicatePESPoints(pes_array, cell, replicate_of=(3, 3) )



######
# DESCRIPTION OF AN INTERFACE OBJECT, FOR DOCUMENTATION
#interface : pymatgen.core.surface.Slab (pymatgen.core.structure.Structure)
        #The interface object of which you want to calculate the PES. It is 
        #needed to extract the lattice parameter used to unfold the energies.
        #type(slab) could be either a Slab object or a Structure object.