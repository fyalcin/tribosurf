#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:56:41 2021

@author: wolloch
"""
from mpinterfaces.transformations import get_aligned_lattices
from pymatgen.core.surface import Slab
from triboflow.firetasks.structure_manipulation import MatchInterface

al111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[2.860165, 0.0, 1.7513459561416446e-16],
   [-1.430082, 2.476975, 1.751345511924877e-16],
   [0.0, 0.0, 32.005944]],
  'a': 2.860165,
  'b': 2.860164274538964,
  'c': 32.005944,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 119.99999682477132,
  'volume': 226.74794103600198},
 'sites': [{'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.0, 0.0, 0.38769624],
   'xyz': [0.0, 0.0, 12.408584146450561],
   'label': 'Al',
   'properties': {'magmom': 0.0}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.333333, 0.666667, 0.46267833],
   'xyz': [-1.0967490001730144e-06, 1.651317492325, 14.80845671999352],
   'label': 'Al',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.0, 0.0, 0.61230376],
   'xyz': [0.0, 0.0, 19.597359853549438],
   'label': 'Al',
   'properties': {'magmom': 0.0}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.666667, 0.333333, 0.53732167],
   'xyz': [1.4300840967489998, 0.8256575076749999, 17.197487280006477],
   'label': 'Al',
   'properties': {'magmom': -0.0}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[2.860165, 0.0, 0.0],
    [-1.430082, 2.476975, 0.0],
    [0.0, 0.0, 32.005944]],
   'a': 2.860165,
   'b': 2.860164274538964,
   'c': 32.005944,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 119.99999682477132,
   'volume': 226.74794103600198},
  'sites': [{'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.0, 0.0, 0.38769624],
    'xyz': [0.0, 0.0, 12.408584146450561],
    'label': 'Al',
    'properties': {'magmom': 0.0}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.333333, 0.666667, 0.46267833],
    'xyz': [-1.0967490001730144e-06, 1.651317492325, 14.80845671999352],
    'label': 'Al',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.0, 0.0, 0.61230376],
    'xyz': [0.0, 0.0, 19.597359853549438],
    'label': 'Al',
    'properties': {'magmom': 0.0}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.666667, 0.333333, 0.53732167],
    'xyz': [1.4300840967489998, 0.8256575076749999, 17.197487280006477],
    'label': 'Al',
    'properties': {'magmom': -0.0}}]},
 'miller_index': (1, 1, 1),
 'shift': 0,
 'scale_factor': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
 'reconstruction': None,
 'energy': None}

cu111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[2.567197, 0.0, 1.571954794415344e-16],
   [-1.2835990000000004, 2.2232580000000004, 1.5719550437332295e-16],
   [0.0, 0.0, 31.288323]],
  'a': 2.567197,
  'b': 2.567197407167007,
  'c': 31.288323,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 120.00000763897559,
  'volume': 178.57939472356944},
 'sites': [{'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.0, -0.0, 0.40101582],
   'xyz': [0.0, 0.0, 12.54711250426986],
   'label': 'Cu',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.333333, 0.666667, 0.46699376],
   'xyz': [-1.6169320002745272e-06, 1.4821727410860004, 14.611451601864479],
   'label': 'Cu',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.0, -0.0, 0.59898418],
   'xyz': [0.0, 0.0, 18.741210495730137],
   'label': 'Cu',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.666667, 0.333333, 0.53300624],
   'xyz': [1.2835996169319999, 0.7410852589140001, 16.676871398135518],
   'label': 'Cu',
   'properties': {'magmom': -0.0}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[2.567197, 0.0, 0.0],
    [-1.283599, 2.223258, 0.0],
    [0.0, 0.0, 31.288323]],
   'a': 2.567197,
   'b': 2.567197407167006,
   'c': 31.288323,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 120.00000763897559,
   'volume': 178.5793947235694},
  'sites': [{'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.0, -0.0, 0.40101582],
    'xyz': [0.0, 0.0, 12.54711250426986],
    'label': 'Cu',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.333333, 0.666667, 0.46699376],
    'xyz': [-1.6169319999784675e-06, 1.482172741086, 14.611451601864479],
    'label': 'Cu',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.0, -0.0, 0.59898418],
    'xyz': [0.0, 0.0, 18.741210495730137],
    'label': 'Cu',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.666667, 0.333333, 0.53300624],
    'xyz': [1.283599616932, 0.741085258914, 16.676871398135518],
    'label': 'Cu',
    'properties': {'magmom': -0.0}}]},
 'miller_index': (1, 1, 1),
 'shift': 0,
 'scale_factor': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
 'reconstruction': None,
 'energy': None}

al111 = Slab.from_dict(al111_dict)
cu111 = Slab.from_dict(cu111_dict)
match_params = {'max_area': 100.0,
                'max_mismatch': 0.1,
                'max_angle_diff': 1.5,
                'r1r2_tol': 0.1
                }

Matcher = MatchInterface(slab_1=al111,
                         slab_2=cu111,
                         strain_factor_1=1,
                         strain_factor_2=1,
                         **match_params)


al111_new, cu111_new = Matcher.get_supercells()

al111_old, cu111_old = get_aligned_lattices(
                al111,
                cu111,
                **match_params)