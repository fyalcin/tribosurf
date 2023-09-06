#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad
from pymatgen.core.structure import Structure

from hitmen_utils.db_tools import VaspDB
from hitmen_utils.vasp_tools import get_custom_vasp_relax_settings
from hitmen_utils.workflows import dynamic_relax_swf

db_file = "/home/mwo/FireWorks/config/db.json"
high_level = "triboflow"
db = VaspDB(db_file=db_file, high_level=high_level)

data = db.find_data(collection="PBE.bulk_data", fltr={"mpid": "mp-30"})

comp_params = data["comp_parameters"]
comp_params["functional"] = "SCAN"
comp_params["use_vdw"] = True

print(comp_params)

prim_bulk = Structure.from_dict(data["structure_equiVol"])

vis = get_custom_vasp_relax_settings(
    prim_bulk, comp_params, "bulk_full_relax"
)

inputs = [[prim_bulk, vis, "Cu_test_SCAN"]]

WF = dynamic_relax_swf(inputs, "test_scan_wf", db_file=db_file)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
