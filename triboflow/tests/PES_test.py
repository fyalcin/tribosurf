#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 14:53:59 2020

@author: wolloch
"""

from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from fireworks import LaunchPad
from triboflow.workflows.subworkflows import CalcPES_SWF
from triboflow.utils.structure_manipulation import SlabFromStructure
from triboflow.utils.database import GetDBJSON, GetSlabFromDB, \
    GetLowEnergyStructure

db_file = GetDBJSON()
functional = "SCAN"

#Test Ag111Ag111 interface
# mpid = "mp-124"
# miller = [1,1,1]
# Slab_dict = GetSlabFromDB(mpid, db_file, miller, functional)
# slab = Slab.from_dict(Slab_dict['relaxed_slab'])
# comp_params = Slab_dict['comp_parameters']

#Test GrapheneGraphene interface
struct, mpid = GetLowEnergyStructure('C', 'mp-1040425')
slab = SlabFromStructure([0,0,1], struct)
comp_params = {'functional': functional,
               'use_vdw': True,
               'use_spin': True,
               'encut': 500,
               'is_metal': True,
               'k_dens': 2000
               }

WF = CalcPES_SWF(top_slab=slab, bottom_slab=slab, comp_parameters=comp_params,
                 file_output=True, output_dir = '/home/fs71332/mwo4',
                 remote_copy = True, server = 'vsc4.vsc.ac.at',  user = 'mwo4', 
                 port = 27)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)