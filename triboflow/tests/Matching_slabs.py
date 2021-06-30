#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:56:41 2021

@author: wolloch
"""
from pymatgen.core.surface import Slab
from mpinterfaces.transformations import get_matching_lattices
from triboflow.utils.database import StructureNavigator

class MatchInterface:
   
    def __init__(self,
                 slab_1,
                 slab_2,
                 strain_factor_1,
                 strain_factor_2,
                 max_area,
                 max_mismatch,
                 max_angle_diff,
                 r1r2_tol):
        self.slab_1 = slab_1.copy()
        self.slab_2 = slab_2.copy()
        self.f1 = strain_factor_1
        self.f2 = strain_factor_2
        self.match_constrains = {'max_area': max_area,
                                 'max_mismatch': max_mismatch,
                                 'max_angle_diff': max_angle_diff,
                                 'r1r2_tol': r1r2_tol
                                 }

    def _get_matching_lattice(self):
        
        uv_opt0, uv_opt1 = get_matching_lattices(self.slab_1,
                                                 self.slab_2,
                                                 **self.match_constrains)
        return uv_opt0, uv_opt1
    

nav = StructureNavigator('auto', True)
al111 = Slab.from_dict(nav.get_slab_from_db('mp-134', 'PBE', [1,1,1])['relaxed_slab'])
cu111 = Slab.from_dict(nav.get_slab_from_db('mp-30', 'PBE', [1,1,1])['relaxed_slab'])
match_params = {'max_area': 100.0,
                'max_mismatch': 0.05,
                'max_angle_diff': 1.5,
                'r1r2_tol': 0.1
                }

Matcher = MatchInterface(slab_1=al111,
                         slab_2=cu111,
                         strain_factor_1=1,
                         strain_factor_2=1,
                         **match_params)
Matcher2 = MatchInterface(slab_1=cu111,
                         slab_2=al111,
                         strain_factor_1=1,
                         strain_factor_2=1,
                         **match_params)

match1 = Matcher._get_matching_lattice()
match2 = Matcher2._get_matching_lattice()

print(match1, match2)
