#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 14:53:59 2020

@author: wolloch
"""

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import calc_pes_swf
from triboflow.utils.structure_manipulation import slab_from_structure
from triboflow.utils.mp_connection import MPConnection
from triboflow.phys.interface_matcher import InterfaceMatcher

functional = "PBE"

mp_connect = MPConnection()
struct, mpid = mp_connect.get_low_energy_structure(
    chem_formula="C", mp_id="mp-1040425"
)
slab = slab_from_structure([0, 0, 1], struct)

IM = InterfaceMatcher(slab, slab, interface_distance=3.4)
interface = IM.get_interface()

comp_params = {
    "functional": functional,
    "use_vdw": True,
    "use_spin": True,
    "encut": 500,
    "is_metal": True,
    "k_dens": 3.0,
}

external_pressure = 1.0

WF = calc_pes_swf(
    interface=interface,
    pressure=external_pressure,
    comp_parameters=comp_params,
    file_output=True,
    output_dir="/fs/home/wolloch/git_test/testdir",
)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
