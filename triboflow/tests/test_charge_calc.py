#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 16:23:42 2021

@author: mwo
"""
from fireworks import LaunchPad
from fireworks import Workflow, Firework, FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from triboflow.utils.database import StructureNavigator

from pymatgen.core.interface import Interface
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.fireworks.core import StaticFW

from triboflow.utils.vasp_tools import get_custom_vasp_static_settings
from triboflow.utils.database import NavigatorMP, Navigator
from triboflow.phys.interface_matcher import InterfaceMatcher
from triboflow.phys.shaper import Shaper


def plot_charge_profile(chgcar, axis=2):
    z=chgcar.zpoints
    x=chgcar.structure.lattice.c*z
    y=chgcar.get_average_along_axis(axis)
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

    
if __name__ == "__main__":
     
    nav_mp = NavigatorMP()
    graphite, _ = nav_mp.get_low_energy_structure('C', mp_id='mp-48')
    nickel, _ = nav_mp.get_low_energy_structure('Ni', mp_id='mp-23')
    
    gr_conv = SpacegroupAnalyzer(graphite).get_conventional_standard_structure()
    ni_conv = SpacegroupAnalyzer(nickel).get_conventional_standard_structure()
    
    sg_params_ni = {'miller': [1,1,1],
                 'slab_thick': 4,
                 'vac_thick': 20,
                 'max_normal_search': 1,
                 'lll_reduce': True,
                 'primitive': True,
                 'tol': 0.1}
    sg_params_gr = {'miller': [0,0,1],
                    'slab_thick': 3,
                    'vac_thick': 20,
                    'max_normal_search': 1,
                    'lll_reduce': True,
                    'primitive': True,
                    'tol': 0.1}
    
    gr_slab_d, _ = Shaper.generate_slabs(gr_conv, sg_params_gr)
    ni_slab_d, _ = Shaper.generate_slabs(ni_conv, sg_params_ni)
    gr_slab = gr_slab_d[(0,0,1)][0]
    ni_slab = ni_slab_d[(1,1,1)][0]
    
    IM = InterfaceMatcher(gr_slab, ni_slab)
    interface = IM.get_interface()
    # FW1 = Firework(tasks=[FT_MakeChargeCalc(structure=slab, 
    #                                        comp_params=comp_params, 
    #                                        calc_name=calc_name)],
    #               name='compute_charge_FW')
    # FW2 = Firework(tasks=[FT_GetCharge(calc_name=calc_name)],
    #               name='get_charge_FW')
    
    # WF = Workflow([FW1, FW2], {FW1: [FW2]}, name='Test_WF_Charge')
    # lpad = LaunchPad.auto_load()
    # lpad.add_wf(WF)
        