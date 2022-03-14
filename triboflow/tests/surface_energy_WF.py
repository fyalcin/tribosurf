#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 11:28:20 2021

@author: yalcin
"""

from fireworks import LaunchPad, Workflow, Firework

from triboflow.firetasks.start_swfs import FT_StartSurfaceEnergySWF
from triboflow.utils.database import Navigator

lpad = LaunchPad.auto_load()

db_file = 'auto'
high_level = 'surfen_test'
override = True
fake = False

nav = Navigator(db_file=db_file, high_level=high_level)
bulk_data = nav.find_many_data('PBE.bulk_data', {})
mat_list = []

for data in bulk_data:
    mpid = data['mpid']
    print(mpid)
    formula = data['formula']
    functional = 'PBE'

    sg_params = {'miller': [(2, 1, 2)],
                 'symmetrize': False,
                 'slab_thick': 12,
                 'vac_thick': 15,
                 'prim': True,
                 'lll_reduce': True,
                 'max_index': 2,
                 'tol': 0.1,
                 'max_normal_search': 'max',
                 'resize': True,
                 'preserve_terminations': True,
                 'min_thick_A': 6}

    sg_filter = {'method': 'all'}

    comp_params_user = {}
    custom_id = None

    SurfaceEnergy = Firework(FT_StartSurfaceEnergySWF(mpid=mpid,
                                                      functional=functional,
                                                      sg_params=sg_params,
                                                      sg_filter=sg_filter,
                                                      db_file=db_file,
                                                      high_level=high_level,
                                                      comp_params_user=comp_params_user,
                                                      custom_id=custom_id),
                             name=f'Start surface energy SWF for {formula}')

    WF = Workflow([SurfaceEnergy])
    formulas = ['WC', 'TiC', 'Fe', 'Pt', 'ZrN', 'TiN', 'GaAs', 'GaN', 'TaN', 'Au', 'Cu', 'ZnCu', 'ZnO', 'CuAu', 'MgO']
    if formula in formulas:
        lpad.add_wf(WF)
