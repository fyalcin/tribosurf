#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from pymatgen.core.structure import Structure
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import converge_kpoints_swf, converge_swf

prim_bulk_dict = {'@module': 'pymatgen.core.structure',
 '@class': 'Structure',
 'charge': None,
 'lattice': {'matrix': [[0.0, 1.98191273760033, 1.98191273760033],
   [1.98191273760033, 0.0, 1.98191273760033],
   [1.98191273760033, 1.98191273760033, 0.0]],
  'a': 2.8028478729543758,
  'b': 2.8028478729543758,
  'c': 2.8028478729543758,
  'alpha': 59.99999999999999,
  'beta': 59.99999999999999,
  'gamma': 59.99999999999999,
  'volume': 15.56981965667947},
 'sites': [{'species': [{'element': 'Pt', 'occu': 1}],
   'abc': [0.0, 0.0, 0.0],
   'xyz': [0.0, 0.0, 0.0],
   'label': 'Pt',
   'properties': {'magmom': 0.000}}]}

struct = Structure.from_dict(prim_bulk_dict)


WF = converge_swf(structure=struct,
                  conv_type='kpoints',
                  flag='NewkpointsConvoTest_SCAN',
                  functional='SCAN',
                  comp_parameters={'is_metal': True,
                                   'use_spin': True,
                                   'functional': 'SCAN'},
                  k_dens_start=2,
                  k_dens_incr=0.1)

WF = converge_swf(structure=struct,
                  conv_type='encut',
                  flag='NewEncutConvoTest_SCAN',
                  functional='SCAN',
                  comp_parameters={'is_metal': True,
                                   'use_spin': True,
                                   'functional': 'SCAN'})


lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
