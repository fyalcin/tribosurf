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
from triboflow.utils.database import GetDBJSON, GetSlabFromDB

db_file = GetDBJSON()
functional = "PBE"
mpid = "mp-124"
miller = [1,1,1]

Slab_dict = GetSlabFromDB(mpid, db_file, miller, functional)
slab = Slab.from_dict(Slab_dict['relaxed_slab'])
comp_params = Slab_dict['comp_parameters']

WF = CalcPES_SWF(top_slab=slab, bottom_slab=slab, comp_parameters=comp_params)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)