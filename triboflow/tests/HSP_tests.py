#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:00:23 2020

@author: wolloch
"""

from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from mpinterfaces.transformations import get_aligned_lattices
from triboflow.helper_functions import GetHighLevelDB


db_file = '/fs/home/wolloch/git_test/config/db.json'

#get an Fe bulk and convert it to a conventional standard structure
db=GetHighLevelDB(db_file)
coll=db.PBE.bulk_data
fe_dict=coll.find_one({'formula': 'Fe'})
fe_bulk=Structure.from_dict(fe_dict.get('structure_equiVol'))
bulk_conv = SpacegroupAnalyzer(fe_bulk).get_conventional_standard_structure()

miller = [1,1,0]
#generate a 110 slab and make a 2x2 supercell of it
SG = SlabGenerator(initial_structure = bulk_conv,
                           miller_index = miller,
                           center_slab = True,
                           primitive = True,
                           max_normal_search=max([abs(l) for l in miller]),
                           lll_reduce = True,
                           min_slab_size = 15,
                           min_vacuum_size = 10)
fe_slab = SG.get_slabs(bonds=None, ftol=0.1, tol=0.1, max_broken_bonds=0,
                            symmetrize=False, repair=False)[0]
fe_2x2_slab=fe_slab.copy()
fe_2x2_slab.make_supercell([2,2,1])

#compute the adsorption sites, add them to the structure and write to file.
adsf = AdsorbateSiteFinder(fe_2x2_slab)
unique_sites = adsf.find_adsorption_sites(distance = 0.5,
                                          symm_reduce = 0.01,
                                          near_reduce = 0.01,
                                          no_obtuse_hollow = True)
fe_2x2_with_ads=fe_2x2_slab.copy()
for s in unique_sites.get('ontop'):
    fe_2x2_with_ads.append('H', s, True)
for s in unique_sites.get('bridge'):
    fe_2x2_with_ads.append('C', s, True)
for s in unique_sites.get('hollow'):
    fe_2x2_with_ads.append('F', s, True)
fe_2x2_with_ads.to(fmt='poscar', filename='fe_2x2_ads.vasp')

#run the same code on a flipped slab
mirror = SymmOp.reflection(normal=[0,0,1], origin=[0, 0, 0])
fe_flipped = fe_2x2_slab.copy()
fe_flipped.apply_operation(mirror)
adsf = AdsorbateSiteFinder(fe_flipped)
unique_sites = adsf.find_adsorption_sites(distance = -0.5,
                                          symm_reduce = 0.01,
                                          near_reduce = 0.01,
                                          no_obtuse_hollow = True)
fe_flipped_ads=fe_flipped.copy()
for s in unique_sites.get('ontop'):
    fe_flipped_ads.append('H', s, True)
for s in unique_sites.get('bridge'):
    fe_flipped_ads.append('C', s, True)
for s in unique_sites.get('hollow'):
    fe_flipped_ads.append('F', s, True)
fe_flipped_ads.to(fmt='poscar', filename='fe_flipped_ads.vasp')

#modify a bottom atom in a 5x5 super cell and flip the structure before computing ads again.
fe_5x5_slab=fe_slab.copy()
fe_5x5_slab.make_supercell([5,5,1])
fe_flipped_ni = fe_5x5_slab.copy()
fe_flipped_ni.replace(12, 'Ni')

fe_flipped_ni.apply_operation(mirror)
adsf = AdsorbateSiteFinder(fe_flipped_ni)
unique_sites = adsf.find_adsorption_sites(distance = -0.5,
                                          symm_reduce = 0.01,
                                          near_reduce = 0.01,
                                          no_obtuse_hollow = True)
fe_flipped_ni_ads=fe_flipped_ni.copy()
for s in unique_sites.get('ontop'):
    fe_flipped_ni_ads.append('H', s, True)
for s in unique_sites.get('bridge'):
    fe_flipped_ni_ads.append('C', s, True)
for s in unique_sites.get('hollow'):
    fe_flipped_ni_ads.append('F', s, True)
fe_flipped_ni_ads.to(fmt='poscar', filename='fe_flipped_ni_ads.vasp')


#finally look at relaxed and unrelaxed slabs and see if there is a difference
db = GetHighLevelDB(db_file)
coll=db.PBE.slab_data
Au_dict = coll.find_one({'formula': 'Au'})
au_relaxed = Structure.from_dict(Au_dict.get('relaxed_slab'))
au_unrelaxed = Structure.from_dict(Au_dict.get('unrelaxed_slab'))
Fe_dict = coll.find_one({'formula': 'Fe', 'miller': [1,1,0]})
fe_unrelaxed = Structure.from_dict(Fe_dict.get('unrelaxed_slab'))
fe_relaxed = Structure.from_dict(Fe_dict.get('relaxed_slab'))

bottom_unrelaxed, top_unrelaxed = get_aligned_lattices(
                fe_unrelaxed,
                au_unrelaxed,
                max_area = 500,
                max_mismatch = 0.1,
                max_angle_diff = 0.1,
                r1r2_tol = 0.1)

bottom_relaxed, top_relaxed = get_aligned_lattices(
                fe_relaxed,
                au_relaxed,
                max_area = 500,
                max_mismatch = 0.1,
                max_angle_diff = 0.1,
                r1r2_tol = 0.1)

adsf = AdsorbateSiteFinder(top_relaxed)
unique_sites_tr = adsf.find_adsorption_sites(distance = 0.5,
                                          symm_reduce = 0.01,
                                          near_reduce = 0.01,
                                          no_obtuse_hollow = True)
adsf = AdsorbateSiteFinder(top_unrelaxed)
unique_sites_tu = adsf.find_adsorption_sites(distance = 0.5,
                                          symm_reduce = 0.01,
                                          near_reduce = 0.01,
                                          no_obtuse_hollow = True)
adsf = AdsorbateSiteFinder(bottom_unrelaxed)
unique_sites_bu = adsf.find_adsorption_sites(distance = 0.5,
                                          symm_reduce = 0.01,
                                          near_reduce = 0.01,
                                          no_obtuse_hollow = True)
adsf = AdsorbateSiteFinder(bottom_relaxed)
unique_sites_br = adsf.find_adsorption_sites(distance = 0.5,
                                          symm_reduce = 0.01,
                                          near_reduce = 0.01,
                                          no_obtuse_hollow = True)