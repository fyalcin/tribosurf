#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:05:37 2021

@author: mwo
"""
from copy import deepcopy
from fireworks import Workflow
from atomate.vasp.fireworks.core import ScanOptimizeFW, OptimizeFW, StaticFW
from atomate.vasp.powerups import add_modify_incar

def scan_relax_with_static(struct, vis, name):
    """Relax a structure with SCAN and do a followup static calculation.
    
    Links a StaticFW to two SanOptimizeFWs. The passed name is given to the
    static calculation, while the PBEsol pre-relaxation will have a
    "_PBEsolPreCalc" postfix, and the SCAN relaxation a "_SCAN_relax" postfix
    to the name.

    Parameters
    ----------
    struct : pymatgen.core.structure.Structure
        structure to relax
    vis : pymatgen.io.vasp.sets.VaspInputSet
        a vasp input set
    name : str
        this will become the task_label in the FireWorks DB. Used to find the
        calculation later.

    Returns
    -------
    list
        list of fireworks
    dict
        dictionary of links between fireworks

    """
    gga_vis = deepcopy(vis)
    static_vis = deepcopy(vis)
    gga_vis.user_incar_settings['ALGO'] = 'Normal'
    vis.user_incar_settings['LWAVE'] = True
    vis_params = {'user_incar_settings': vis.user_incar_settings,
                  'user_kpoints_settings': vis.user_kpoints_settings,
                  'user_potcar_functional': vis.potcar_functional,
                  'vdw': vis.vdw}
    fw_1 = ScanOptimizeFW(structure=struct, name=name+'_PBEsolPreCalc',
                          vasp_input_set=gga_vis,
                          vasp_input_set_params={'vdw': vis.vdw})
    fw_2 = ScanOptimizeFW(vasp_input_set_params=vis_params,
                          parents = fw_1, prev_calc_loc=True,
                          name=name+'_SCAN_relax')
    static_vis.user_incar_settings['EDIFF'] = 0.5e-06
    static_vis.user_incar_settings['ISMEAR'] = -5
    static_vis.user_incar_settings['ALGO'] = 'Normal'
    static_vis.user_incar_settings['SIGMA'] = 0.05
    static_vis.user_incar_settings['ISTART'] = 1
    static_vis.user_incar_settings['ICHARG'] = 0
    static_vis.user_incar_settings['NELMDL'] = 0
    static_vis_params = {'user_incar_settings': static_vis.user_incar_settings,
                         'user_kpoints_settings': static_vis.user_kpoints_settings,
                         'user_potcar_functional': static_vis.potcar_functional,
                         'user_potcar_settings': static_vis.user_potcar_settings,
                         'vdw': static_vis.vdw}
    fw_3 = StaticFW(vasp_input_set_params=static_vis_params,
                              parents = fw_2, prev_calc_loc=True,
                              additional_files_from_prev_calc=['WAVECAR'],
                              name=name)
    return [fw_1, fw_2, fw_3], {fw_1: fw_2, fw_2: fw_3}

def pbe_relax_with_static(struct, vis, name):
    """Relax a structure with PBE and do a followup static calculation.
    
    Links a StaticFW to an OptimizeFW. The passed name is given to the static
    calculation, while the relaxation will have a "_PBE_relax" postfix to the
    name.

    Parameters
    ----------
    struct : pymatgen.core.structure.Structure
        structure to relax
    vis : pymatgen.io.vasp.sets.VaspInputSet
        a vasp input set
    name : str
        this will become the task_label in the FireWorks DB. Used to find the
        calculation later.

    Returns
    -------
    list
        list of fireworks
    dict
        dictionary of links between fireworks

    """
    vis.user_incar_settings['LWAVE'] = True
    static_vis = deepcopy(vis)
    static_vis.user_incar_settings['ISMEAR'] = -5
    static_vis.user_incar_settings['EDIFF'] = 1e-07
    static_vis.user_incar_settings['SIGMA'] = 0.05
    static_vis.user_incar_settings['ISTART'] = 1
    static_vis.user_incar_settings['ICHARG'] = 0
    static_vis.user_incar_settings['NELMDL'] = 0
    vis_params = {'user_incar_settings': static_vis.user_incar_settings,
                  'user_kpoints_settings': vis.user_kpoints_settings,
                  'user_potcar_functional': vis.potcar_functional,
                  'user_potcar_settings': vis.user_potcar_settings,
                  'vdw': vis.vdw}

    fw_1 = OptimizeFW(structure=struct,
                    name = name+'_PBE_relax',
                    vasp_input_set=vis)
    fw_2 = StaticFW(vasp_input_set_params=vis_params,
                    parents = fw_1, prev_calc_loc=True,
                    additional_files_from_prev_calc=['WAVECAR'],
                    name=name)
    return [fw_1, fw_2], {fw_1: fw_2}

def scan_relax(struct, vis, name):
    """Relax a structure with SCAN.
    
    Links two SanOptimizeFWs. The passed name is given to the second one,
    which is a SCAN relaxation, while the PBEsol pre-relaxation will have a
    "_PBEsolPreCalc" postfix.

    Parameters
    ----------
    struct : pymatgen.core.structure.Structure
        structure to relax
    vis : pymatgen.io.vasp.sets.VaspInputSet
        a vasp input set
    name : str
        this will become the task_label in the FireWorks DB. Used to find the
        calculation later.

    Returns
    -------
    list
        list of fireworks
    dict
        dictionary of links between fireworks

    """
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
    return [fw_1, fw_2], {fw_1: fw_2}

def pbe_relax(struct, vis, name):
    """Relax a structure with PBE.
    
    Use an OptimizeFW to relax a strucutre. The name will be assigned to the
    task_label of the calculation.

    Parameters
    ----------
    struct : pymatgen.core.structure.Structure
        structure to relax
    vis : pymatgen.io.vasp.sets.VaspInputSet
        a vasp input set
    name : str
        this will become the task_label in the FireWorks DB. Used to find the
        calculation later.

    Returns
    -------
    list
        list of fireworks
    dict
        dictionary of links between fireworks

    """
    fw = OptimizeFW(structure=struct, name=name, vasp_input_set=vis)
    return [fw], {}

def dynamic_relax_swf(inputs_list,
                      wf_name = 'Dynamicaly generated Workflow',
                      add_static = False):
    """Generate a workflow with relaxations for PBE or SCAN functionals.
    
    A base workflow that uses an input list of lists with structures, vasp-
    input-sets and tags (to find the calculation again in the Fireworks
    database) to create a workflow. Many structures can be optimized at once.
    If the 'SCAN' (or '2rSCAN') is set as METAGGA value for any of the vasp-
    input-sets, two ScanOptimizeFW are run in succession for each structure,
    if not, a single OptimizeFW is run. A add_modify_incar powerup from
    atomate is applied to the final workflow. Static calculations can be
    added in the end if the 'add_static' flag is set to True. This will result
    in the input names pointing to the static follow up run. Those static runs
    will start from the WAVECAR of the previous relaxation and provide more
    accurate total energies.
    

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
    add_static : bool, optional
        Selects if a static follow-up calculation should be run at the relaxed
        positions starting from the CHGCAR of the relaxation run. The default
        is False

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A workfow intended to relax a structure or many structures.

    """
    fw_list = []
    fw_links = {}
    SCAN_list = ['scan', 'rscan', 'r2scan','Scan', 'Rrscan', 'R2scan',
                 'SCAN', 'RSCAN', 'R2SCAN']
    for i in range(len(inputs_list)):
        struct = inputs_list[i][0]
        vis = inputs_list[i][1]
        name = inputs_list[i][2]
        if add_static:
            if vis.incar.get('METAGGA') in SCAN_list:
                fws, links = scan_relax_with_static(struct, vis, name)
                fw_list.extend(fws)
                fw_links.update(links)
            else:
                fws, links = pbe_relax_with_static(struct, vis, name)
                fw_list.extend(fws)
                fw_links.update(links)
        else:
            if vis.incar.get('METAGGA') in SCAN_list:
                fws, links = scan_relax(struct, vis, name)
                fw_list.extend(fws)
                fw_links.update(links)
            else:
                fws, links = pbe_relax(struct, vis, name)
                fw_list.extend(fws)
                fw_links.update(links)
    WF = Workflow(fireworks=fw_list,
                  links_dict=fw_links,
                  name=wf_name)
    return add_modify_incar(WF)