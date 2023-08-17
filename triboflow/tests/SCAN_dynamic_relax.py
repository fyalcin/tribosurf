#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from pymatgen.core.structure import Structure
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.workflows.base import dynamic_relax_swf
from triboflow.utils.database import StructureNavigator
from hitmen_utils.vasp_tools import get_custom_vasp_relax_settings

db_file = "/home/mwo/FireWorks/config/db.json"

nav_structure = StructureNavigator(db_file=db_file, high_level="triboflow")
data = nav_structure.get_bulk_from_db(mp_id="mp-30", functional="PBE")

comp_params = data["comp_parameters"]
comp_params["functional"] = "SCAN"
comp_params["use_vdw"] = True

print(comp_params)

prim_bulk = Structure.from_dict(data["structure_equiVol"])

vis = get_custom_vasp_relax_settings(prim_bulk, comp_params, "bulk_full_relax")

inputs = [[prim_bulk, vis, "Cu_test_SCAN"]]

WF = dynamic_relax_swf(inputs, "test_scan_wf")

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
