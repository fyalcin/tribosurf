#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from pymatgen.core.structure import Structure
from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from atomate.vasp.fireworks.core import StaticFW
from triboflow.helper_functions import GetBulkFromDB


db_file = '/home/mwo/FireWorks/config/db.json'

data = GetBulkFromDB("mp-81", db_file, 'PBE')

struct = Structure.from_dict(data['structure_equiVol'])
compdat

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
rapidfire(lpad)
