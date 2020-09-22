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
# UTILITY FUNCTIONS
# =============================================================================


def plot_slab_atoms_hs(material, hs, lattice, element='element', mill='mill'):    
 
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from mpl_toolkits.mplot3d import Axes3D    
    
    hs = np.array(hs)
    a, b = lattice
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_slab(material, ax, repeat=5, draw_unit_cell=True, adsorption_sites=True)
    plt.plot(hs[:,0], hs[:,1], 'o')
    ax.set(xlim = (-np.sign(a)*0.5, hetero.lattice.a+np.sign(a)*0.5), 
           ylim = (-np.sign(b)*0.5, hetero.lattice.b+np.sign(b)*0.5))
    plt.savefig(element1+mill1_str+'_hs.pdf', dpi=300)


def plot_slab_hs_red(hs, lattice, out='hs.pdf', c='k', tot=False):
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from mpl_toolkits.mplot3d import Axes3D
        
    a, b = lattice
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rect = patches.Rectangle((0, 0), a, b, linewidth=1, edgecolor='k',
                             facecolor='none')
    ax.add_patch(rect)
    if tot:
        sub = np.array(hs[0])
        top = np.array(hs[1])
        plt.plot(sub[:,0], sub[:,1], 'o', color='b')
        plt.plot(top[:,0], top[:,1], 'o', color='r')
    else:
        hs = np.array(hs)
        plt.plot(hs[:,0], hs[:,1], 'o', color=c)
    ax.set_aspect('equal')
    plt.savefig(out, dpi=300)
    plt.show()
    

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
                                     max_normal_search = None,
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

# During a copy the x axis is oriented along 100, if the basal plane of the
# cell is not orthogonal we have problem when running get_interface

substrate = substrate.copy()
coating   = coating.copy()


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

def ExtractHSPoints(material, near_reduce=0.01):
    
    # COLLECT UNIQUE HS POINTS
    AdsMat = AdsorbateSiteFinder(material)
    unique_sites = AdsMat.find_adsorption_sites( distance = 0, 
                                       symm_reduce = 1e-06, 
                                       near_reduce = near_reduce, 
                                       no_obtuse_hollow = True )
    
    #COLLECT ALL THE HS POINTS
    replica_sites = AdsMat.find_adsorption_sites( distance = 0, 
                                       symm_reduce = 0, 
                                       near_reduce = near_reduce, 
                                       no_obtuse_hollow = True )
    
    # Extract and identify each unique HS point of materials and its replica
    HS = GetUniqueHS(unique_sites)             
    HS_ALL = GetReplicaHS(AdsMat, replica_sites, HS)

    
    return HS, HS_ALL

def GetUniqueHS(unique_sites):
    
    HS = {} 
    for k in ['ontop', 'bridge', 'hollow']:
        if unique_sites[k] != []:
            n=1
            for data in unique_sites[k]:
                HS[k+'_'+str(n)] = data
                n += 1
    HS['all'] = unique_sites['all']
    
    return HS
    

def GetReplicaHS(AdsMat, replica_sites, hs):
    
    # b is the object obtained from the AdsorbateSiteFinder function
    # b_sites is the corresponding list of sites, replicated inside whole cell    
    # Create an empty dictionary with the list of unique HS points as keys
    HS_ALL = {}
    
    for k in hs.keys():
        HS_ALL[k] = []
    
    # Recognize of which unique HS point each sites is the replica
    for key in HS_ALL.keys():
        for site in replica_sites['all']:        
            
            HSToEvaluate = [hs[key].copy()]
            HSToEvaluate.append(np.array(site))
            
            HSEvaluated = AdsMat.symm_reduce(HSToEvaluate)
          
            if len(HSEvaluated) == 1:
                HS_ALL[key].append(site)
    
    HS_ALL['all'] = replica_sites['all']

    return HS_ALL


hs_sub, hs_sub_all = ExtractHSPoints(substrate, 0.01)
hs_coat, hs_coat_all = ExtractHSPoints(coating, 0.01)


# =============================================================================
# PLOT AND PRINT INFORMATION
# =============================================================================


# Extract information

print(element1+mill1_str+' replica: '+str(hetero.lattice.a/slab1.lattice.a)+
      'X' + str(hetero.lattice.b/slab1.lattice.b))
print(element2+mill2_str+' replica: '+str(hetero.lattice.a/slab2.lattice.a)+
      'X' + str(hetero.lattice.b/slab2.lattice.b))


# Plot and save the aligned cells displaying the HS points
# WARNING: plot_slab of pymatgen call adsorption_sites with symm_reduce!=0

plot_slab_atoms_hs(substrate, hs_sub['all'], (hetero.lattice.a, hetero.lattice.b), 
             element1, mill1_str)

plot_slab_atoms_hs(coating, hs_coat['all'], (hetero.lattice.a, hetero.lattice.b), 
             element2, mill2_str)


# Plot the HS in the cell -- Replicated HS (!!!)

plot_slab_hs_red(hs_sub_all['all'], (hetero.lattice.a, hetero.lattice.b), 
                 element1+mill1_str+'_hs_all.pdf', c='b')
plot_slab_hs_red(hs_coat_all['all'], (hetero.lattice.a, hetero.lattice.b), 
                 element2+mill2_str+'hs_all.pdf', c='r')

plot_slab_hs_red((hs_sub['all'], hs_coat['all']), (hetero.lattice.a, hetero.lattice.b), 
                 'superposition_hs.pdf', c='r', tot=True)


# =============================================================================
# CALCULATE COMBINATIONS OF HS POINTS - TEST
# =============================================================================

# Only working with squared lattices, at the moment
def combine_hs(hetero, hs_sub, hs_coat): 
    
    alat1 = hetero.lattice.a
    alat2 = hetero.lattice.b
    
    hs_inter = []
    sites = []
        
    ksub  = list(hs_sub.keys());  ksub.remove('all')
    kcoat = list(hs_coat.keys()); kcoat.remove('all')
            
    for k1 in ksub:
        for k2 in kcoat:
                
            for i, d1 in enumerate(np.matrix(list(hs_sub[k1]))):
                for j, d2 in enumerate(np.matrix(list(hs_coat[k2]))):
                    
                    #print(k1, d1, d2)
                    shift = np.array(d1 - d2)
                    #print(shift)
                    
                    if shift[0,0] < 0:
                        shift[0,0] += alat1
                    elif shift[0,0] >= alat1:
                        shift[0,0] -= alat1
                    if shift[0,1] < 0:
                        shift[0,1] += alat2
                    elif shift[0,1] >= alat2:
                        shift[0,1] -= alat2
                        
                    hs_inter.append([shift[0,0], shift[0,1]])
                    sites.append('sub-'+k1+'#'+str(i)+'_VS_''coat-'+k2+'#'+str(j))
    return np.array(hs_inter), sites

hs_inter, sites = combine_hs(hetero, hs_sub, hs_coat)
plot_slab_hs_red(hs_inter, (hetero.lattice.a, hetero.lattice.b), 'unfolded_hs.pdf',
                 c='k')

hs_inter_1, sites_1 = combine_hs(hetero, hs_sub, hs_coat_all)
plot_slab_hs_red(hs_inter_1, (hetero.lattice.a, hetero.lattice.b), 'unfolded1_hs.pdf',
                 c='k')

hs_inter_2, sites_2 = combine_hs(hetero, hs_sub_all, hs_coat)
plot_slab_hs_red(hs_inter_2, (hetero.lattice.a, hetero.lattice.b), 'unfolded2_hs.pdf',
                 c='k')


# =============================================================================
# CALCULATE SHIFTS TO BE DONE FOR THE PES
# =============================================================================




# =============================================================================
# TEST TO RECOGNIZE THE HS POINTS
# =============================================================================

# # 1 version ###
# bridge1 = []
# bridge2 = []

# for site in b_sites['bridge']:
#     site_list = b.symm_reduce([hs_sub['bridge_1'], site])
#     if len(site_list) == 1:
#         bridge1.append(site)
#     else:
#         bridge2.append(site)
# ###      

# # 2nd partial version

# dict_final = {}
# for k in hs_subs.keys():
#     dict_final[k] = 0
# dict_final.pop('all')

# for symm_points in ['bridge', 'hollow', 'ontop']:
#    for site in b_sites[symm_points]:
#        LIST = []
       
#        for key in dict_final.keys():
#            site_list = b.symm_reduce([hs_sub[key], site])
#            if len(site_list) == 1:
#                bridge1.append(site)
#            else:
#                bridge2.append(site)


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
