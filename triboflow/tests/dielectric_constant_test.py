#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:29:36 2021

@author: mwo
"""

from uuid import uuid4

from fireworks import LaunchPad

from pymatgen.core.structure import Structure

from triboflow.utils.database import StructureNavigator
from triboflow.workflows.dielectric import dielectric_constant_swf


mpid = 'mp-1265'
functional = 'PBE'

nav = StructureNavigator('auto', True)
struct_dict = nav.get_bulk_from_db(mpid, functional)
struct = Structure.from_dict(struct_dict['structure_equiVol'])
comp_params = struct_dict['comp_parameters']

flag = 'test_dielectric_WF_'+str(uuid4())

WF = dielectric_constant_swf(struct, mpid, flag, comp_params,
                             update_bulk=False, update_slabs=True)

lpad = LaunchPad().auto_load()
lpad.add_wf(WF)