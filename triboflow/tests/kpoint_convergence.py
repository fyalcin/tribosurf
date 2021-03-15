#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from pymatgen.core.structure import Structure
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import converge_kpoints_swf
from triboflow.utils.database import GetBulkFromDB


db_file = '/home/mwo/FireWorks/config/db.json'

#data = GetBulkFromDB("mp-134", db_file, 'PBE') #Al
data = GetBulkFromDB("mp-81", db_file, 'PBE') #Au
#data = GetBulkFromDB("mp-66", db_file, 'PBE') #diamond C

struct = Structure.from_dict(data['structure_fromMP'])
comp_parameters = data['comp_parameters']

comp_parameters['encut'] = 500

WF = converge_kpoints_swf(struct, comp_parameters, {}, data['mpid'],
                          comp_parameters['functional'], k_dens_start=500,
                          k_dens_incr=50)

#mod_WF = add_modify_incar(WF, modify_incar_params={'incar_update': uis})

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
rapidfire(lpad)
