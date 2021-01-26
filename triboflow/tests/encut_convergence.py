#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from pymatgen.core.structure import Structure
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import ConvergeEncut_SWF
from triboflow.helper_functions import GetBulkFromDB


db_file = '/home/mwo/FireWorks/config/db.json'

data = GetBulkFromDB("mp-134", db_file, 'PBE')

struct = Structure.from_dict(data['structure_fromMP'])
comp_parameters = data['comp_parameters']

comp_parameters['encut'] = 500

WF = ConvergeEncut_SWF(struct, comp_parameters, {}, data['mpid'],
                         comp_parameters['functional'])


lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
