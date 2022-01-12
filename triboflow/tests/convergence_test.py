#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from pymatgen.core.structure import Structure
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import converge_swf

prim_bulk_dict = {'@module': 'pymatgen.core.structure',
 '@class': 'Structure',
 'charge': None,
 'lattice': {'matrix': [[-1.420025838669831,
    1.4200258386698312,
    1.4200258386698312],
   [1.420025838669831, -1.4200258386698312, 1.4200258386698312],
   [1.4200258386698315, 1.4200258386698312, -1.420025838669831]],
  'a': 2.4595569006367533,
  'b': 2.4595569006367533,
  'c': 2.4595569006367533,
  'alpha': 109.47122063449069,
  'beta': 109.47122063449069,
  'gamma': 109.47122063449069,
  'volume': 11.453777224502772},
 'sites': [{'species': [{'element': 'Fe', 'occu': 1}],
   'abc': [0.0, 0.0, 0.0],
   'xyz': [0.0, 0.0, 0.0],
   'label': 'Fe',
   'properties': {'magmom': 2.214}}]}

struct = Structure.from_dict(prim_bulk_dict)


# WF = converge_swf(structure=struct,
#                   conv_type='kpoints',
#                   flag='NewkspacingConvoTest_SCAN',
#                   functional='SCAN',
#                   comp_parameters={'is_metal': True,
#                                    'use_spin': True,
#                                    'functional': 'SCAN'},
#                   k_dens_start=2,
#                   k_dens_incr=0.1)

WF = converge_swf(structure=struct,
                  conv_type='encut',
                  flag='SCAN_encut_test_Fe_100',
                  functional='SCAN',
                  encut_incr= 100,
                  comp_parameters={'is_metal': True,
                                    'use_spin': True,
                                    'functional': 'SCAN',
                                    'k_dens': 4.6})


lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
