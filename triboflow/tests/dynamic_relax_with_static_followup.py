#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 11:17:06 2021

@author: wolloch
"""

from fireworks import LaunchPad

from pymatgen.core.surface import Slab
from pymatgen.core.structure import Structure
from triboflow.workflows.base import dynamic_relax_swf
from triboflow.utils.database import StructureNavigator
from triboflow.utils.vasp_tools import get_custom_vasp_static_settings, get_custom_vasp_relax_settings

nav = StructureNavigator('auto', True)

mpids = ['mp-134', 'mp-13', 'mp-66']
input_list = []
for mpid in mpids:
    mat_dict = nav.get_bulk_from_db(functional='PBE', mp_id=mpid)
    struct = Structure.from_dict(mat_dict['primitive_structure'])
    comp_params = mat_dict['comp_parameters']
    comp_params['functional'] = 'SCAN'
    #comp_params['use_vdw'] = True
    vis = get_custom_vasp_relax_settings(struct, comp_params, 'bulk_full_relax')
    input_list.append([struct, vis, f'{struct.composition.reduced_formula}_test_relax_SCAN'])

WF = dynamic_relax_swf(inputs_list=input_list, wf_name='test_relax_wf',
                       add_static=True)
lpad = LaunchPad.auto_load()
lpad.add_wf(WF)