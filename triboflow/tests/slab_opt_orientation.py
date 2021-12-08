#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 11:28:20 2021

@author: yalcin
"""

from fireworks import Firework, Workflow
from fireworks import LaunchPad
from triboflow.firetasks.surfen import FT_SlabOptOrientation

mpid = 'mp-149'
functional = 'PBE'
db_file = 'auto'
high_level = True
comp_params = {}

max_index = 3

override = True

sg_params = {'miller': [(1, 1, 0)],
             'symmetrize': False,
             'slab_thick': 15,
             'vac_thick': 15,
             'prim': True,
             'lll_reduce': True,
             'minimize_bv': True,
             'tol': 0.1,
             'mns': 1}

sg_filter = {'method': 'all'}

bvs_method = 'all'
bvs_param = 1
fake = True

lpad = LaunchPad.auto_load()

# for lll in [True, False]:
#     for prim in [True, False]:
#         for mns in [1, 'max']:
#             for sym in [True, False]:
#                 sg_params = {'symmetrize': sym,
#                              'slab_thick': 15,
#                              'vac_thick': 15,
#                              'prim': prim,
#                              'lll_reduce': lll,
#                              'minimize_bv': True,
#                              'tol': 0.1,
#                              'mns': mns}
#
#                 FT = FT_SlabOptOrientation(mpid=mpid,
#                                            functional=functional,
#                                            db_file=db_file,
#                                            max_index=max_index,
#                                            high_level=high_level,
#                                            comp_params=comp_params,
#                                            bvs_method=bvs_method,
#                                            sg_params=sg_params,
#                                            bvs_param=bvs_param,
#                                            override=override,
#                                            fake=fake)
#
#                 FW = Firework(FT)
#                 WF1 = Workflow.from_Firework(FW)
#
#                 lpad.add_wf(WF1)

FT = FT_SlabOptOrientation(mpid=mpid,
                           functional=functional,
                           sg_params=sg_params,
                           sg_filter=sg_filter,
                           db_file=db_file,
                           high_level=high_level,
                           comp_params=comp_params,
                           override=override,
                           fake=fake)

FW = Firework(FT)
WF1 = Workflow.from_Firework(FW)

lpad.add_wf(WF1)
# rapidfire(lpad)
