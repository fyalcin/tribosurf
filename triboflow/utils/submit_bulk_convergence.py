#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 16:12:32 2021

@author: wolloch
"""

from fireworks import LaunchPad
from fireworks import Workflow, Firework
from triboflow.fireworks.init_fws import InitWF
from triboflow.firetasks.structure_manipulation import FT_StartPreRelax
from triboflow.firetasks.init_check import material_from_mp
from triboflow.firetasks.convergence import FT_StartConvo
from triboflow.utils.utils import load_homoatomic_materials

computational_params = {'functional': 'PBE',
                   'volume_tolerance': 0.001,
                   'BM_tolerance': 0.01,
                   'use_vdw': False,
                   'surfene_thr': 0.01,
                   'vacuum': 15}

def get_bulk_convergence_wf(material, comp_params):
    struct, mp_id = material_from_mp(material)
    functional = comp_params.get('functional', 'PBE')
    
    WF = []

    Initialize = InitWF.checkinp_bulk_convo(material=material,
                                            computational=comp_params)
    WF.append(Initialize)

    PreRelaxation = Firework(FT_StartPreRelax(mp_id=mp_id,
                                              functional=functional),
                                 name='Start pre-relaxation for {}'
                                .format(material['formula']))
    WF.append(PreRelaxation)
    
    ConvergeEncut = Firework(FT_StartConvo(conv_type='encut',
                                           mp_id=mp_id,
                                           functional=functional,
                                           ),
                                name='Start encut convergence for {}'
                                      .format(material['formula']))
    WF.append(ConvergeEncut)
    
    ConvergeKpoints = Firework(FT_StartConvo(conv_type='kpoints',
                                           mp_id=mp_id,
                                           functional=functional,
                                           ),
                                name='Start kpoints convergence for {}'
                                      .format(material['formula']))
    WF.append(ConvergeKpoints)
    
    Dependencies = {Initialize: [PreRelaxation],
                    PreRelaxation: [ConvergeEncut],
                    ConvergeEncut: [ConvergeKpoints]}

    WF_Name = 'ConvergeBulk ' + material['formula'] +' '+ mp_id +' '+ functional

    WF = Workflow(WF, Dependencies, name=WF_Name)
    
    return WF

def submit_multiple_wfs(workflow_list):    
    lpad = LaunchPad.auto_load()
    for wf in workflow_list:
        lpad.add_wf(wf)



if __name__ == "__main__":
    materials_dict = load_homoatomic_materials()
    materials_list =[]

    for formula in ['Al', 'C', 'Si', 'Ge', 'Cu', 'Ag',
                    'Au', 'Ni', 'Fe', 'Ti', 'Co']:
        mpid = materials_dict[formula]['mpids'][materials_dict[formula]['default']]
        materials_list.append({'formula': formula, 'mpid': mpid})
        workflow_list = []
    
    for mat in materials_list:
        workflow_list.append(get_bulk_convergence_wf(mat, computational_params))
    submit_multiple_wfs(workflow_list)
