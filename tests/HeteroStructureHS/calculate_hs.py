#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" SCRIPT PROTOTYPE FOR THE PROCEDURE TO BE FOLLOWED IN THE WORKFLOW

Created on Tue Jul 14 16:00:00 2020
@author: glosi000
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, Slab
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.transformations.standard_transformations import RotationTransformation

from mpinterfaces.transformations import get_aligned_lattices, get_interface
#from mpinterfaces.transformations import get_aligned_lattices, generate_all_configs


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

area = 300         # Max possible area for the hetero interface
mismatch = 0.02    # Max mismatch between the lattice parameters of the cells
angle_diff = 0.1   # Max mismatch between orientation of lattice parameters
r1r2 = 0.2         # Tolerance

n1 = thickness1    # Number of layers of slab1 in the interface
n2 = thickness2    # Number of layers of slab2 in the interface
sep = 2            # Interface separation in Angstrom


# =============================================================================
# GENERATE THE STRUCTURE AND THE SLABS
# =============================================================================


# Read structure from files element.POSCAR
# WARNING: CONVENTIONAL UNIT CELLS ARE NEEDED TO GENERATE THE SLABS

mat1 = Structure.from_file(element1+'.POSCAR')
mat1 = SpacegroupAnalyzer(mat1).get_conventional_standard_structure()

mat2 = Structure.from_file(element2+'.POSCAR')
mat2 = SpacegroupAnalyzer(mat2).get_conventional_standard_structure()


folder=element2+mill2_str+'_on_'+element1+mill1_str
try:
    os.mkdir(folder)
except:
    pass
os.chdir(folder)


# Generate the slabs

mill1 = [ int(i) for i in list(mill1_str) ]
mill2 = [ int(i) for i in list(mill2_str) ]


mat1_SlabGenerator = SlabGenerator(  mat1, mill1, 
                                     min_slab_size     = thickness1, 
                                     min_vacuum_size   = vacuum,
                                     lll_reduce        = lll_reduce, 
                                     center_slab       = centered, 
                                     in_unit_planes    = unit_plane, 
                                     primitive         = primitive, 
                                     max_normal_search = 5,
                                     reorient_lattice  = True      )

slab1 = mat1_SlabGenerator.get_slab()
slab1.to('poscar', element1 + mill1_str + '.POSCAR')


if layered == False :    # Slab is generated only if mat2 is not a 2D material
    
    mat2_SlabGenerator = SlabGenerator(  mat2, mill2, 
                                         min_slab_size     = thickness2, 
                                         min_vacuum_size   = vacuum,
                                         lll_reduce        = lll_reduce, 
                                         center_slab       = centered, 
                                         in_unit_planes    = unit_plane, 
                                         primitive         = primitive, 
                                         max_normal_search = 5,
                                         reorient_lattice  = True      )
    
    slab2 = mat2_SlabGenerator.get_slab()
    slab2.to(fmt = 'poscar', filename = element2 + mill2_str + '.POSCAR')

else:
    
    slab2 = mat2


# =============================================================================
# CREATE THE HETERO INTERFACE
# =============================================================================


substrate, coating = get_aligned_lattices( slab1, slab2,
                                           max_area       = area,
                                           max_mismatch   = mismatch,
                                           max_angle_diff = angle_diff,
                                           r1r2_tol       = r1r2 )

# Convert Structures to Slab python type

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


substrate.to(fmt = 'poscar', filename = element1 + mill1_str + '_align.POSCAR')
coating.to(fmt = 'poscar', filename = element2 + mill2_str + '_align.POSCAR')


hetero = get_interface( substrate, coating, 
                        nlayers_substrate      = n1,
                        nlayers_2d             = n2, 
                        separation             = sep )

Poscar(hetero).write_file(filename='POSCAR')    

# OLD VERSION
# hetero = generate_all_configs(coating, substrate, n2, n1,sep)
# for i, iface in enumerate(hetero):
#     Poscar(iface).write_file(filename='hetero.POSCAR_{}'.format(i))
    

# =============================================================================
# CALCULATE HS POINTS
# =============================================================================


a = AdsorbateSiteFinder(substrate)
subs_hs = np.stack(a.find_adsorption_sites( distance=0, 
                                            symm_reduce=0, near_reduce=0.01, 
                                            no_obtuse_hollow=True ) ['all'] )

# MPInterfaces rotates the slab when building the interface
b = AdsorbateSiteFinder(coating)
coat_hs = np.stack(b.find_adsorption_sites( distance=0, 
                                            symm_reduce=0, near_reduce=0.01, 
                                            no_obtuse_hollow=True ) ['all'] )


# Extract information

print(element1+mill1_str+' replica: '+str(hetero.lattice.a/slab1.lattice.a)+
      'X' + str(hetero.lattice.b/slab1.lattice.b))
print(element2+mill2_str+' replica: '+str(hetero.lattice.a/slab2.lattice.a)+
      'X' + str(hetero.lattice.b/slab2.lattice.b))


# Plot and save the aligned cells displaying the HS points
# WARNING: plot_slab of pymatgen call adsorption_sites with symm_reduce!=0

fig = plt.figure()
ax = fig.add_subplot(111)
plot_slab(substrate, ax, repeat=5, draw_unit_cell=True, adsorption_sites=True)
plt.plot(subs_hs[:,0], subs_hs[:,1], 'o')
ax.set(xlim = (-np.sign(hetero.lattice.a)*0.5, hetero.lattice.a+np.sign(hetero.lattice.a)*0.5), 
       ylim = (-np.sign(hetero.lattice.b)*0.5, hetero.lattice.b+np.sign(hetero.lattice.b)*0.5))
plt.savefig(element1+mill1_str+'_hs.pdf', dpi=300)

fig = plt.figure()
ax = fig.add_subplot(111)
plot_slab(coating, ax, repeat=5, draw_unit_cell=True, adsorption_sites=True)
plt.plot(coat_hs[:,0], coat_hs[:,1], 'o')
ax_max = max(hetero.lattice.a, hetero.lattice.b)
ax.set(xlim = (-np.sign(hetero.lattice.a)*0.5, 
               hetero.lattice.a + np.sign(hetero.lattice.a)*0.5), 
       ylim = (-np.sign(hetero.lattice.b)*0.5,
               hetero.lattice.b + np.sign(hetero.lattice.b)*0.5))
plt.savefig(element2+mill2_str+'_hs.pdf', dpi=300)


# Plot the HS in the cell

fig = plt.figure()
ax = fig.add_subplot(111)
rect = patches.Rectangle((0, 0), hetero.lattice.a, hetero.lattice.b,
                         linewidth=1, edgecolor='k',facecolor='none')
ax.add_patch(rect)
plt.plot(subs_hs[:,0], subs_hs[:,1], 'o', color='b')
ax.set_aspect('equal')
plt.savefig(element1+mill1_str+'_hs_all.pdf', dpi=300)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
rect = patches.Rectangle((0, 0), hetero.lattice.a, hetero.lattice.b,
                         linewidth=1, edgecolor='k',facecolor='none')
ax.add_patch(rect)
plt.plot(coat_hs[:,0], coat_hs[:,1], 'o', color='r')
ax.set_aspect('equal')
plt.savefig(element2+mill2_str+'_hs_all.pdf', dpi=300)

fig = plt.figure()
ax = fig.add_subplot(111)
rect = patches.Rectangle((0, 0), hetero.lattice.a, hetero.lattice.b,
                         linewidth=1, edgecolor='k',facecolor='none')
ax.add_patch(rect)
plt.plot(subs_hs[:,0], subs_hs[:,1], 'o', color='b')
plt.plot(coat_hs[:,0], coat_hs[:,1], 'o', color='r')
ax.set_aspect('equal')
plt.savefig('superposition_hs_all.pdf', dpi=300)


# =============================================================================
# CALCULATE COMBINATIONS OF HS POINTS
# =============================================================================


all_sites = []



# =============================================================================
# NOTES / OLD CODE
# =============================================================================


# NOT WORKING: HS points of the slabs (MPINTERFACES WITHOUT MODIFICATION)

# a = SpacegroupAnalyzer(substrate, symprec=0.1, angle_tolerance=5.0)
# a = a.get_conventional_standard_structure()
# a2 = AdsorbateSiteFinder(a)
# subs_hs = np.stack(a2.find_adsorption_sites( distance=0, 
#                                              symm_reduce=0, near_reduce=0, 
#                                              no_obtuse_hollow=False ) ['all'] )


# b = SpacegroupAnalyzer(coating, symprec=0.1, angle_tolerance=5.0)
# b = b.get_conventional_standard_structure()
# b2 = AdsorbateSiteFinder(b)
# coat_hs = np.stack(b2.find_adsorption_sites( distance=0, 
#                                              symm_reduce=0, near_reduce=0, 
#                                              no_obtuse_hollow=False ) ['all'] )