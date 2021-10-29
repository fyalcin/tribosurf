#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 17:00:07 2021

@author: mwo
"""

from fireworks import FWAction, FiretaskBase, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.fireworks.core import StaticFW

from triboflow.utils.vasp_tools import get_custom_vasp_static_settings
from triboflow.utils.database import Navigator


def plot_charge_profile(chgcar, axis=2):
    z=chgcar.zpoints
    x=chgcar.structure.lattice.c*z
    y=chgcar.get_average_along_axis(2)
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(x, y)
    fig.savefig(f'charge_profile_{chgcar.structure.formula}.png')
    
    

@explicit_serialize
class FT_MakeChargeCalc(FiretaskBase):
    
    required_params = ['structure', 'comp_params', 'calc_name']
    optional_params = ['db_file', 'high_level_db']
    
    def run_task(self, fw_spec):
        
        struct = self.get('structure')
        comp_params = self.get('comp_params')
        label = self.get('calc_name')
        
        vis = get_custom_vasp_static_settings(struct, comp_params,
                                              'slab_from_scratch')
                       
        FW = StaticFW(structure=struct,
                      vasp_input_set=vis,
                      name=label,
                      vasptodb_kwargs = {'store_volumetric_data': ['chgcar']})
        WF = add_modify_incar(Workflow.from_Firework(FW))
        return FWAction(detours=WF)
    
@explicit_serialize
class FT_GetCharge(FiretaskBase):
    
    required_params = ['calc_name']
    optional_params = ['db_file', 'high_level_db']
    
    def run_task(self, fw_spec):
        label = self.get('calc_name')
        nav = Navigator()
        chgcar = nav.get_chgcar_from_label(label)
        plot_charge_profile(chgcar)
        
        
        