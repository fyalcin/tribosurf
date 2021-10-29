#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 16:23:42 2021

@author: mwo
"""
from fireworks import LaunchPad
from fireworks import Workflow, Firework


from triboflow.utils.database import StructureNavigator
from triboflow.firetasks.charge_test import FT_GetCharge, FT_MakeChargeCalc

from pymatgen.core.surface import Slab
from pymatgen.core.structure import Structure

nav = StructureNavigator('auto', True)
slab_dict = nav.get_slab_from_db('mp-134', 'PBE', [1,0,0])
# slab = Slab.from_dict(slab_dict.get('relaxed_slab'))
slab = Structure.from_dict(slab_dict['thickness']['data_12']['output']['structure'])
comp_params = slab_dict.get('comp_parameters')
calc_name = f'interface_{slab.formula}'

    
if __name__ == "__main__":
     
    FW1 = Firework(tasks=[FT_MakeChargeCalc(structure=slab, 
                                           comp_params=comp_params, 
                                           calc_name=calc_name)],
                  name='compute_charge_FW')
    FW2 = Firework(tasks=[FT_GetCharge(calc_name=calc_name)],
                  name='get_charge_FW')
    
    WF = Workflow([FW1, FW2], {FW1: [FW2]}, name='Test_WF_Charge')
    lpad = LaunchPad.auto_load()
    lpad.add_wf(WF)
        