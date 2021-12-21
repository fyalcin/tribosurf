#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 11:28:20 2021

@author: yalcin
"""

from fireworks import LaunchPad

from triboflow.workflows.subworkflows import surface_energy_swf

mpid = 'mp-149'
functional = 'PBE'
db_file = 'auto'
high_level = True
comp_params = {}
override = True
fake = False

sg_params = {'miller': [(1, 1, 1)],
             'symmetrize': False,
             'slab_thick': 8,
             'vac_thick': 15,
             'prim': True,
             'lll_reduce': True,
             'minimize_bv': True,
             # 'max_index': 2,
             'tol': 0.1,
             'max_normal_search': 'max'}

sg_filter = {'method': 'all'}

lpad = LaunchPad.auto_load()

WF = surface_energy_swf(mpid=mpid,
                        functional=functional,
                        sg_params=sg_params,
                        sg_filter=sg_filter,
                        db_file=db_file,
                        high_level=high_level,
                        comp_params_user=comp_params)

lpad.add_wf(WF)
