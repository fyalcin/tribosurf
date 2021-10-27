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
from triboflow.firetasks.check_inputs import FT_CopyCompParamsToSlab
from triboflow.firetasks.init_check import material_from_mp
from triboflow.firetasks.convergence import FT_StartConvo
from triboflow.utils.utils import load_homoatomic_materials
from triboflow.firetasks.run_slabs_wfs import FT_SlabOptThick

computational_params = {'functional': 'PBE',
                   'volume_tolerance': 0.001,
                   'BM_tolerance': 0.01,
                   'use_vdw': False,
                   'surfene_thr': 0.01,
                   'vacuum': 15}

def get_slab_convergence_wf(material, comp_params):
    struct, mp_id = material_from_mp(material)
    formula = material['formula']
    functional = comp_params.get('functional', 'PBE')
    miller=material.get('miller')
    FT_CopyCompParamsToSlab
    WF = []

    Initialize = InitWF.checkinp_slab_convo(material=material,
                                            computational=comp_params)
    WF.append(Initialize)

    PreRelaxation = Firework(FT_StartPreRelax(mp_id=mp_id,
                                              functional=functional),
                                 name='Start pre-relaxation for {}'
                                .format(formula))
    WF.append(PreRelaxation)
    
    ConvergeEncut = Firework(FT_StartConvo(conv_type='encut',
                                           mp_id=mp_id,
                                           functional=functional,
                                           ),
                                name='Start encut convergence for {}'
                                      .format(formula))
    WF.append(ConvergeEncut)
    
    ConvergeKpoints = Firework(FT_StartConvo(conv_type='kpoints',
                                           mp_id=mp_id,
                                           functional=functional,
                                           ),
                                name='Start kpoints convergence for {}'
                                      .format(formula))
    WF.append(ConvergeKpoints)
    
    TransferParams = Firework(FT_CopyCompParamsToSlab(mp_id=mp_id,
                                                      miller=miller,
                                                      functional=functional),
                              name='Copy converged params')
    WF.append(TransferParams)
    
    MakeSlabs = Firework(FT_SlabOptThick(mp_id=mp_id,
                                         miller=miller,
                                         functional=functional),
                         name=f'Slab thickness optimization for {formula}')
    WF.append(MakeSlabs)
    
    Dependencies = {Initialize: [PreRelaxation],
                    PreRelaxation: [ConvergeEncut],
                    ConvergeEncut: [ConvergeKpoints],
                    ConvergeKpoints: [TransferParams],
                    TransferParams: [MakeSlabs]}

    WF_Name = f'ConvergeSlab {formula} {mp_id} {miller} {functional}'

    WF = Workflow(WF, Dependencies, name=WF_Name)
    
    return WF

def submit_multiple_wfs(workflow_list):    
    lpad = LaunchPad.auto_load()
    for wf in workflow_list:
        lpad.add_wf(wf)



if __name__ == "__main__":
    materials_dict = load_homoatomic_materials()
    materials_list =[]
    miller_list = ['100', '110', '111']

    # for formula in ['Al', 'C', 'Si', 'Ge', 'Cu', 'Ag',
    #                 'Au', 'Ni', 'Fe', 'Ti', 'Co']:
    for formula in ['Al', 'C']:
        mpid = materials_dict[formula]['mpids'][materials_dict[formula]['default']]
        mat_dict = {'formula': formula, 'mpid': mpid}
        for m in miller_list:
            mat_dict['miller'] = m
            mat_dict['thick_min'] = 4
            mat_dict['thick_max'] = 12
            mat_dict['thick_incr'] = 1
            if formula in ['C', 'Si', 'Ge'] and m == '111':
                mat_dict['thick_incr'] = 2
            materials_list.append(mat_dict.copy())
        
    workflow_list = []
    
    for mat in materials_list:
        workflow_list.append(get_slab_convergence_wf(mat, computational_params))
    submit_multiple_wfs(workflow_list)
