#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from fireworks import LaunchPad
from pymatgen.core.structure import Structure

from hitmen_utils.db_tools import VaspDB
# from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import converge_swf

db = VaspDB("auto", True)
bulk_dict = db.find_data(collection="PBE.bulk_data", fltr={"mpid": "mp-30"})
struct = Structure.from_dict(bulk_dict["structure_fromMP"])
comp_parameters = bulk_dict["comp_parameters"]
mpid = bulk_dict["mpid"]

WF = converge_swf(
    structure=struct,
    conv_type="encut",
    flag=f"encut_test for {mpid}",
    comp_parameters=comp_parameters,
)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
