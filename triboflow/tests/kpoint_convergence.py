#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.io.vasp.inputs import Kpoints
from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.powerups import add_modify_incar
from triboflow.helper_functions import (
        GetBulkFromDB, GetCustomVaspStaticSettings)


db_file = '/home/mwo/FireWorks/config/db.json'

data = GetBulkFromDB("mp-81", db_file, 'PBE')

struct = Structure.from_dict(data['structure_equiVol'])
comp_parameters = data['comp_parameters']

vis, uis, vdw = GetCustomVaspStaticSettings(struct, comp_parameters,
                                            'bulk_from_scratch')
kpoints = Kpoints.automatic_density(struct, 100)
print(kpoints.as_dict())
vis = MPStaticSet(struct, user_incar_settings=uis,
                  user_kpoints_settings=kpoints)
FW = StaticFW(structure=struct, vasp_input_set=vis, name='unique_label')

WF = Workflow.from_Firework(FW, name='test WF')

#mod_WF = add_modify_incar(WF, modify_incar_params={'incar_update': uis})

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
rapidfire(lpad)
