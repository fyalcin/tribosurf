#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Created on Wed Apr 15 17:04:35 2020

@author: mwo
"""

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from CommonWorkflows import Relax_SWF
from HelperFunctions import GetLowEnergyStructure, GetGapFromMP

# =============================================================================
# THIS IS A SIMPLE TEST AS A COMPLETE WORKFLOW:

# struct , mp_id = GetLowEnergyStructure('NiO')

# prev_spec = {'some_data': {'a': 1, 'b':2},
#              'KPOINTS': 'previous_KPTS',
#              'relax_loc': './RELAX_THIS_THING'}

# bandgap = GetGapFromMP(mp_id)
# if bandgap > 0.1:
#     is_metal = False
# else:
#     is_metal = True

# parameters = {'use_vdw': False,
#               'functional': 'PBE',
#               #'k_distance': 25,
#               #'encut': 350,
#               'is_metal': is_metal}

# Relaxed_structure_loc = ['structures',
#                          'relaxed_structures',
#                          struct.composition.reduced_formula]
# wf = Relax_SWF(struct, parameters, 'bulk_pos_relax', Relaxed_structure_loc,
#                prev_spec)
# #fw = OptimizeFW(struct, override_default_vasp_params=uis)
# #fw2 = RelaxFW(name='print spec')

# lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
# #lpad.add_wf(Workflow([fw, fw2], {fw: [fw2]}))
# lpad.add_wf(wf)

# rapidfire(lpad)

# =============================================================================
# THIS IS A TEST WHERE A THE SUBWORKFLOW IS CALLED FROM ANOTHER WORKFLOW:
# =============================================================================

from atomate.vasp.powerups import clean_up_files
from CommonFiretasks import FT_PrintSpec, FT_FetchStructureFromFormula
from CommonFireworks import StartDetourWF_FW
from CommonFireworks import CheckInputsFW



#Create the input dictionary

inputs = {'material_1': {'formula': 'Al',
                          'miller': '111',
                          'min_vacuum': 25,
                          'min_thickness': 5
                          },
          'material_2': {'formula': 'Fe',
                          'mp_id': 'mp-13',
                          'miller': '110',
                          'min_vacuum': 25,
                          'min_thickness': 6
                          },
          'computational_params':{'functional': 'PBE',
                                  'use_vdw': 'True'},
          'interface_params':{'interface_distance': 1.5,
                              'max_area': 500,
                              'r1r2_tol': 0.05
                              }
          }

# create the Workflow List
FWs = []


Initialize = Firework(FT_PrintSpec(), spec=inputs,
                                name= 'Initialize Workflow')
FWs.append(Initialize)

Check_Inputs = CheckInputsFW(mat1loc=['material_1'], mat2loc=['material_2'],
                             compparamloc=['computational_params'],
                             interparamloc=['interface_params'],
                             name='Check Inputs')
FWs.append(Check_Inputs)

Start_M1 = Firework(FT_FetchStructureFromFormula(
                        materials_dict_loc=['material_1']),
                        name='Get '+inputs['material_1']['formula']+' bulk')
FWs.append(Start_M1)

Start_M2 = Firework(FT_FetchStructureFromFormula(
                        materials_dict_loc=['material_2']),
                        name='Get '+inputs['material_2']['formula']+' bulk')
FWs.append(Start_M2)

bulk_loc = ['material_1', inputs['material_1']['formula']+'_fromMP']
out_loc = ['material_1', inputs['material_1']['formula']+'_relaxed']
Relax_M1 = StartDetourWF_FW('Relax_SWF',
                            name='Relax '+inputs['material_1']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_1'],
                            relax_type='bulk_full_relax',
                            out_loc=out_loc,
                            to_pass=['material_1'])
FWs.append(Relax_M1)

bulk_loc = ['material_2', inputs['material_2']['formula']+'_fromMP']
out_loc = ['material_2', inputs['material_2']['formula']+'_relaxed']
Relax_M2 = StartDetourWF_FW('None',
                            name='Relax '+inputs['material_2']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_2'],
                            relax_type='bulk_full_relax',
                            out_loc=out_loc,
                            to_pass=['material_2'])
FWs.append(Relax_M2)

End = Firework(FT_PrintSpec(), name='End of the workflow')
FWs.append(End)

Dependencies = {Initialize: [Check_Inputs],
                Check_Inputs: [Start_M1, Start_M2],
                Start_M1: [Relax_M1],
                Start_M2: [Relax_M2],
                Relax_M1: [End],
                Relax_M2: [End]}

wf = Workflow(FWs, Dependencies, name='Test Relax Subworkflow')
mod_wf = clean_up_files(wf, ('WAVECAR*', 'CHGCAR*'))

# finally, instatiate the LaunchPad and add the workflow to it
lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(mod_wf)
rapidfire(lpad)