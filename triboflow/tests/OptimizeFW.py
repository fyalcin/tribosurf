#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire#
from atomate.vasp.fireworks.core import OptimizeFW
from triboflow.fireworks.common import CheckInputsFW
from triboflow.helper_functions import GetLowEnergyStructure
from triboflow.helper_functions import GetCustomVaspRelaxSettings, GetBulkFromDB

db_file='/home/mwo/FireWorks/config/db.json'

data = GetBulkFromDB('mp-81', db_file, 'PBE')
miller = [1, 1, 1]
min_thickness = 10
min_vacuum = 20

comp_params = data['comp_parameters']
comp_params['k_dens'] = 1000

prim_bulk = Structure.from_dict(data['structure_equiVol'])
conv_bulk = SpacegroupAnalyzer(prim_bulk).get_conventional_standard_structure()
        
SG = SlabGenerator(initial_structure = conv_bulk,
                   miller_index = miller,
                   min_slab_size = min_thickness,
                   min_vacuum_size = min_vacuum)
slab = SG.get_slabs(bonds=None, ftol=0.1, tol=0.1, max_broken_bonds=0,
                    symmetrize=False, repair=False)[0]

vis = GetCustomVaspRelaxSettings(slab, comp_params, 'slab_pos_relax')

Optimize_FW = OptimizeFW(slab, name='Test_relaxation',
                         vasp_input_set = vis,
                         half_kpts_first_relax = True)
WF = Workflow([Optimize_FW], name='OptimizeFW WF')

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
rapidfire(lpad)
