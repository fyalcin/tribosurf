#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:05:37 2021

@author: mwo
"""
from copy import deepcopy
from fireworks import Workflow
from atomate.vasp.fireworks.core import ScanOptimizeFW, OptimizeFW
from atomate.vasp.powerups import add_modify_incar

def dynamic_relax_swf(inputs_list,
                      wf_name = 'Dynamicaly generated Workflow'):
    """Generate a workflow with relaxations for PBE or SCAN functionals.
    
    A base workflow that uses an input list of lists with structures, vasp-
    input-sets and tags (to find the calculation again in the Fireworks
    database) to create a workflow. Many structures can be optimized at once.
    If the 'SCAN' (or '2rSCAN') is set as METAGGA value for any of the vasp-
    input-sets, two ScanOptimizeFW are run in succession for each structure,
    if not, a single OptimizeFW is run. A add_modify_incar powerup from
    atomate is applied to the final workflow.
    

    Parameters
    ----------
    inputs_list : list of lists
        List containing further list of the type:
            [[structure_1, vasp_input_set_1, name_1],
             [structure_2, vasp_input_set_2, name_2],
             ....
             [structure_n, vasp_input_set_n, name_n]]
    wf_name : str, optional
        Name for the created workflow. The default is 'Dynamicaly generated
        Workflow'.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A workfow intended to relax a structure or many structures.

    """
    fw_list = []
    fw_links = {}
    for i in range(len(inputs_list)):
        struct = inputs_list[i][0]
        vis = inputs_list[i][1]
        name = inputs_list[i][2]
        if vis.incar.get('METAGGA') in ['Scan', 'R2scan']:
            gga_vis = deepcopy(vis)
            gga_vis.user_incar_settings['ALGO'] = 'Normal'
            vis_params = {'user_incar_settings': vis.user_incar_settings,
                          'user_kpoints_settings': vis.user_kpoints_settings,
                          'user_potcar_functional': vis.potcar_functional,
                          'vdw': vis.vdw}
            fw_1 = ScanOptimizeFW(structure=struct, name=name+'_PBEsolPreCalc',
                                  vasp_input_set=gga_vis,
                                  vasp_input_set_params={'vdw': vis.vdw})
            fw_2 = ScanOptimizeFW(vasp_input_set_params=vis_params,
                                      parents = fw_1, prev_calc_loc=True,
                                      name=name)
        
            fw_list.extend([fw_1, fw_2])
            fw_links[fw_1] = fw_2
        else:
            fw = OptimizeFW(structure=struct,
                            name = name,
                            vasp_input_set=vis)
            fw_list.append(fw)
            
    WF = Workflow(fireworks=fw_list,
                  links_dict=fw_links,
                  name=wf_name)
    return add_modify_incar(WF)