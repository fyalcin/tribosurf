#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Created on Wed Apr 15 17:04:35 2020

@author: mwo
"""

import numpy as np
import pprint
from uuid import uuid4
from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.workflows.base.bulk_modulus import get_wf_bulk_modulus
from CommonFiretasks import FT_PrintSpec, FT_UpdateBMLists, \
    FT_EnergyCutoffConvo
from HelperFunctions import GetLowEnergyStructure, GetCustomVaspStaticSettings

# =============================================================================
# THIS IS A SIMPLE TEST AS A COMPLETE WORKFLOW:

# struct , mp_id = GetLowEnergyStructure('Al')

# deforms = []
# for i in np.arange(0.9, 1.1, 0.05):
#     dm=np.eye(3)*i
#     deforms.append(dm)

# vis, uis, vdw = GetCustomVaspStaticSettings(struct, {}, 'bulk_from_scratch')

# pprint.pprint(uis)

# wf = get_wf_bulk_modulus(struct, deforms, vasp_input_set=None,
#                          vasp_cmd=VASP_CMD,
#                          db_file=DB_FILE, user_kpoints_settings=None,
#                          eos='birch_murnaghan', tag=None,#'test_BM',
#                          user_incar_settings=uis)

# FW = Firework([FT_UpdateBMLists(formula=struct.composition.reduced_formula),
#               FT_PrintSpec()], name='Print BM List')

# wf.append_wf(Workflow.from_Firework(FW), wf.leaf_fw_ids)

# lpad = LaunchPad.auto_load()
# lpad.add_wf(wf)


# rapidfire(lpad)


# =============================================================================
# Try out the SWF in a propper setting.
# =============================================================================

struct , mp_id = GetLowEnergyStructure('Au')
tag = "BM group: {}".format(str(uuid4()))

FT1 = FT_EnergyCutoffConvo(structure = struct, out_loc = ['here'],
                          comp_params = {'is_metal': True},
                          tag = tag,
                          encut_incr = 25,
                          encut_start = 200)
FT2 = FT_PrintSpec()

FW = Firework([FT1, FT2], name='Converge Encut Firework')

WF = Workflow.from_Firework(FW, name='Encut Convergence Workflow')

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
rapidfire(lpad)

# from atomate.vasp.powerups import clean_up_files
# from CommonFiretasks import FT_PrintSpec, FT_FetchStructureFromFormula
# from CommonFireworks import CheckInputsFW, StartDetourWF_FW



# #Create the input dictionary

# inputs = {'material_1': {'formula': 'FeO',
#                           'miller': '111',
#                           'min_vacuum': 25,
#                           'min_thickness': 5
#                           },
#           'material_2': {'formula': 'Al',
#                           #'mp_id': 'mp-13',
#                           'miller': '110',
#                           'min_vacuum': 25,
#                           'min_thickness': 6
#                           },
#           'computational_params':{'functional': 'PBE',
#                                   'use_vdw': 'False'},
#           'interface_params':{'interface_distance': 1.5,
#                               'max_area': 500,
#                               'r1r2_tol': 0.05
#                               }
#           }

# # create the Workflow List
# FWs = []


# Initialize = Firework(FT_PrintSpec(), spec=inputs,
#                                 name= 'Initialize Workflow')
# FWs.append(Initialize)

# Check_Inputs = CheckInputsFW(mat1loc=['material_1'], mat2loc=['material_2'],
#                              compparamloc=['computational_params'],
#                              interparamloc=['interface_params'],
#                              name='Check Inputs')
# FWs.append(Check_Inputs)

# Start_M1 = Firework(FT_FetchStructureFromFormula(
#                         materials_dict_loc=['material_1']),
#                         name='Get '+inputs['material_1']['formula']+' bulk')
# FWs.append(Start_M1)

# Start_M2 = Firework(FT_FetchStructureFromFormula(
#                         materials_dict_loc=['material_2']),
#                         name='Get '+inputs['material_2']['formula']+' bulk')
# FWs.append(Start_M2)

# bulk_loc = ['material_1', inputs['material_1']['formula']+'_fromMP']
# out_loc = ['material_1']
# KPTSConvo_M1 = StartDetourWF_FW('Converge_Kpoints_SWF',
#                             name='Converge Kpoints '+
#                             inputs['material_1']['formula'],
#                             structure_loc=bulk_loc,
#                             comp_parameters_loc=['material_1'],
#                             out_loc=out_loc,
#                             to_pass=['material_1'])
# FWs.append(KPTSConvo_M1)

# bulk_loc = ['material_2', inputs['material_2']['formula']+'_fromMP']
# out_loc = ['material_2']
# KPTSConvo_M2 = StartDetourWF_FW('Converge_Kpoints_SWF',
#                             name='Converge Kpoints '+
#                             inputs['material_2']['formula'],
#                             structure_loc=bulk_loc,
#                             comp_parameters_loc=['material_2'],
#                             out_loc=out_loc,
#                             to_pass=['material_2'])
# FWs.append(KPTSConvo_M2)

# End = Firework(FT_PrintSpec(), name='End of the workflow')
# FWs.append(End)

# Dependencies = {Initialize: [Check_Inputs],
#                 Check_Inputs: [Start_M1, Start_M2],
#                 Start_M1: [KPTSConvo_M1],
#                 Start_M2: [KPTSConvo_M2],
#                 KPTSConvo_M1: [End],
#                 KPTSConvo_M2: [End]}

# wf = Workflow(FWs, Dependencies, name='Test Kpoints Convergence Subworkflow')
# mod_wf = clean_up_files(wf, ('WAVECAR*', 'CHGCAR*'))

# # finally, instatiate the LaunchPad and add the workflow to it
# lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
# lpad.add_wf(mod_wf)
# #rapidfire(lpad)