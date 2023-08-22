#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 14:53:59 2020

@author: wolloch
"""

from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import calc_pes_swf
from triboflow.utils.structure_manipulation import slab_from_structure
from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.mp_connection import MPConnection
from triboflow.phys.interface_matcher import InterfaceMatcher

# nav = Navigator()
# db_file = nav.path
functional = "PBE"

# TODO: TO REFACTOR
# Test Ag111Ag111 interface
# mpid = "mp-124"
# miller = [1,1,1]
# nav_structure = StructureNavigator(
#     db_file=db_file,
#     high_level='triboflow')
# slab_dict = nav_structure.get_slab_from_db(
#     mp_id=mpid,
#     functional=functional,
#     miller=miller)
# slab = Slab.from_dict(slab_dict['relaxed_slab'])
# comp_params = slab_dict['comp_parameters']

# Test Au111_C001 interface
# int_dict = GetInterfaceFromDB('Au111_C001_mp-66_mp-81', db_file, functional)
# top_slab = Slab.from_dict(int_dict['top_aligned'])
# bottom_slab = Slab.from_dict(int_dict['bottom_aligned'])
# comp_params = int_dict['comp_parameters']


# WF = CalcPES_SWF(top_slab=top_slab, bottom_slab=bottom_slab, top_mpid='mp-66',
#                  bottom_mpid='mp-81',comp_parameters=comp_params,
#                  file_output=True, output_dir = '/home/fs71332/mwo4',
#                  remote_copy = True, server = 'vsc4.vsc.ac.at',  user = 'mwo4',
#                  port = 27)

# lpad = LaunchPad.auto_load()
# lpad.add_wf(WF)
# Test GrapheneGraphene interface
mp_connect = MPConnection()
struct, mpid = mp_connect.get_low_energy_structure(chem_formula="C", mp_id="mp-1040425")
slab = slab_from_structure([0, 0, 1], struct)

IM = InterfaceMatcher(slab, slab, interface_distance=3.4)
interface = IM.get_interface()

# inter_dict = interface.as_dict()
# bottom_dict = slab.as_dict()
# top_dict = slab.as_dict()

# nav_high = Navigator('auto', high_level=True)
# nav_high.update_data(
#     collection=functional + '.interface_data',
#     fltr={'name': 'Navigator_test_Gr-Gr'},
#     new_values={'$set': {'unrelaxed_structure': inter_dict,
#                          'bottom_aligned': bottom_dict,
#                          'top_aligned': top_dict}})

comp_params = {
    "functional": functional,
    "use_vdw": True,
    "use_spin": True,
    "encut": 500,
    "is_metal": True,
    "k_dens": 3.0,
}

external_pressure = 1.0

WF = calc_pes_swf(interface=interface, pressure=external_pressure, comp_parameters=comp_params, file_output=True,
                  output_dir="/fs/home/wolloch/git_test/testdir")

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
