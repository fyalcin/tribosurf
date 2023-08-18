#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from pymatgen.core.structure import Structure
from fireworks import LaunchPad

# from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import converge_swf
from triboflow.utils.database import StructureNavigator


# db_file = "/home/mwo/FireWorks/config/db.json"

# nav_structure = StructureNavigator(db_file=db_file, high_level="triboflow")
# data = nav_structure.get_bulk_from_db(mp_id="mp-30", functional="PBE")

# struct = Structure.from_dict(data["structure_fromMP"])
# comp_parameters = data["comp_parameters"]


WF = converge_swf(
    struct, comp_parameters, {}, data["mpid"], comp_parameters["functional"]
)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
