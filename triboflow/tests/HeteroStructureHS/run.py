#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Test the modules written for HS points and PES/shear strength
Created on Mon Oct  5 11:27:30 2020

@author: gl
"""

import os
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import SlabGenerator, Slab
from mpinterfaces.transformations import get_aligned_lattices, get_interface
from HS_functions import *
from utility_functions import *


# =============================================================================
# INPUT
# =============================================================================


element1 = 'Au'
mill1_str = '100'

element2 = 'V'
mill2_str = '110'


# =============================================================================
# CONTROL PANEL
# =============================================================================


# Paramteres for the SLABS
# WARNING: CENTERED SHOULD BE SET TO TRUE OR MPINTERFACES GIVES PROBLEMS

layered = False    # Wether element2 is a 2D Material
thickness1 = 6     # Number of atomic layers / Min Angstrom thicknessin slab1
thickness2 = 6     # Number of atomic layers / Min Angstrom thicknessin slab2
vacuum = 25        # Number of atomic layers for vacuum / Angstrom
lll_reduce = True  # FORCE THE C VECTOR TO BE ORTHOGONAL TO THE XY PLANE
primitive = False  # Reduce the slab cell to its primitive to reduce atoms
unit_plane = True  # Set thick_slab, vacuum in units of planes or Angstrom
centered = True    # Wether to center the slab at th middle of the cell


# Parameters for the INTERFACE

area = 500         # Max possible area for the hetero interface
mismatch = 0.05    # Max mismatch between the lattice parameters of the cells
angle_diff = 0.1   # Max mismatch between orientation of lattice parameters
r1r2 = 0.2         # Tolerance

n1 = thickness1    # Number of layers of slab1 in the interface
n2 = thickness2    # Number of layers of slab2 in the interface
sep = 2            # Interface separation in Angstrom


# =============================================================================
# MAIN
# =============================================================================


if __name__ == '__main__':
    
    ###########################################################################
    ########################### Load the structures ###########################
    ###########################################################################
    
    os.chdir('./Structures')
    
    mat1 = Structure.from_file(element1+'.POSCAR')
    mat1 = SpacegroupAnalyzer(mat1).get_conventional_standard_structure()

    mat2 = Structure.from_file(element2+'.POSCAR')
    mat2 = SpacegroupAnalyzer(mat2).get_conventional_standard_structure()
    
    mill1 = [ int(i) for i in list(mill1_str) ]
    mill2 = [ int(i) for i in list(mill2_str) ]
    
    
    ###########################################################################
    ###########################  Create the slabs  ############################
    ###########################################################################
    
    
    mat1_SlabGenerator = SlabGenerator(  mat1, mill1, 
                                     min_slab_size     = thickness1, 
                                     min_vacuum_size   = vacuum,
                                     lll_reduce        = lll_reduce, 
                                     center_slab       = centered, 
                                     in_unit_planes    = unit_plane, 
                                     primitive         = primitive, 
                                     max_normal_search = None,
                                     reorient_lattice  = True      )
    
    mat2_SlabGenerator = SlabGenerator(  mat2, mill2, 
                                         min_slab_size     = thickness2, 
                                         min_vacuum_size   = vacuum,
                                         lll_reduce        = lll_reduce, 
                                         center_slab       = centered, 
                                         in_unit_planes    = unit_plane, 
                                         primitive         = primitive, 
                                         max_normal_search = 5,
                                         reorient_lattice  = True      )
    
    slab1 = mat1_SlabGenerator.get_slab()
    slab2 = mat2_SlabGenerator.get_slab()
    
    
    ###########################################################################
    #########################  Create the interface  ##########################
    ###########################################################################
    
    
    substrate, coating = get_aligned_lattices( slab1, slab2,
                                              max_area       = area,
                                              max_mismatch   = mismatch,
                                              max_angle_diff = angle_diff,
                                              r1r2_tol       = r1r2 )
    
    substrate = Slab(substrate.lattice,
                      substrate.species_and_occu,
                      substrate.frac_coords,
                      mill1,
                      Structure.from_sites(substrate, to_unit_cell=True),
                      shift=0,
                      scale_factor=np.eye(3, dtype=np.int),
                      site_properties=substrate.site_properties)

    coating = Slab(coating.lattice,
                    coating.species_and_occu,
                    coating.frac_coords,
                    mill2,
                    Structure.from_sites(coating, to_unit_cell=True),
                    shift=0,
                    scale_factor=np.eye(3, dtype=np.int),
                    site_properties=coating.site_properties)
    
    substrate = substrate.copy()
    coating   = coating.copy()
    
    
    hetero = get_interface( substrate, coating, 
                            nlayers_substrate      = n1,
                            nlayers_2d             = n2, 
                            separation             = sep )
    
    Poscar(hetero).write_file(filename='POSCAR')
    
    
    ###########################################################################
    ##########################  Get the HS points  ############################
    ###########################################################################
    
    
    # Extract the HS points for the slabs
    hs_sub, hs_sub_all = GetSlabHS(substrate)
    hs_coat, hs_coat_all = GetSlabHS(coating)
    
    # Calculate the HS points for the hetero interface
    hs = GetInterfaceHS(hs_sub, hs_coat)
    hs_all = GetInterfaceHS(hs_sub_all, hs_coat_all)
    
    hs = ApplyPbcToHS(hetero, hs)
    hs_all = ApplyPbcToHS(hetero, hs_all)
    
    # Remove the z-coordinates 
    hs_sub = RemoveZCoords(hs_sub)
    hs_sub_all = RemoveZCoords(hs_sub_all)
    hs_coat = RemoveZCoords(hs_coat)
    hs_coat_all = RemoveZCoords(hs_coat_all)
    hs = RemoveZCoords(hs)
    hs_all = RemoveZCoords(hs_all)
    
    # Plot the HS points
    Plot_SlabHS(substrate,  hs_sub,  to_fig=None)
    Plot_SlabHS(  coating,  hs_coat, to_fig=None)
    Plot_SlabHS(substrate,       hs, to_fig=None)