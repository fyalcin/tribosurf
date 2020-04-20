#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:37:34 2020

@author: mwo
"""

import os
from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from CommonFiretasks import *
from CommonFireworks import *

indir = os.getcwd()


#Create the input dictionary

inputs = {'material_1': {'formula': 'Cu',
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
                                  'use_vdw': 'False',
                                  'use_spin': 'True'
                                  },
          'interface_params':{'interface_distance': 1.5,
                              'max_area': 500,
                              'r1r2_tol': 0.05
                              }
          }

# create the Workflow List
WF = []

Initialize = Firework(FT_PrintSpec(), spec=inputs,
                                name= 'Initialize Workflow')
WF.append(Initialize)

Check_Inputs = CheckInputsFW(mat1name='material_1', mat2name='material_2',
                             compparaname='computational_params',
                             interparaname='interface_params',
                             name='Check Inputs')
WF.append(Check_Inputs)

Start_M1 = Firework(FT_FetchStructureFromFormula(
                        formula_loc=['material_1','formula'],
                        structure_name=inputs['material_1']['formula']+'bulk'),
                        name='Get '+inputs['material_1']['formula']+' bulk')
WF.append(Start_M1)

Start_M2 = Firework(FT_FetchStructureFromFormula(
                        formula_loc=['material_2','formula'],
                        structure_name=inputs['material_2']['formula']+'bulk'),
                        name='Get '+inputs['material_2']['formula']+' bulk')
WF.append(Start_M2)

Converge_M1 = ConvergeParametersFW(name='Converge Parameters for'+' '+
                                        inputs['material_1']['formula'])
WF.append(Converge_M1)

Converge_M2 = ConvergeParametersFW(name='Converge Parameters for'+' '+
                                        inputs['material_2']['formula'])
WF.append(Converge_M2)

Final_Params = FixParametersFW(name='Select computational parameters')
WF.append(Final_Params)

Relax_M1 = RelaxFW(name='Relax '+inputs['material_1']['formula']+' bulk')
WF.append(Relax_M1)

Relax_M2 = RelaxFW(name='Relax '+inputs['material_2']['formula']+' bulk')
WF.append(Relax_M2)

Make_Slab_M1 = Firework(FT_MakeSlabFromStructure(
                        bulk_name=inputs['material_1']['formula']+'bulk',
                        dict_name='material_1'),
                name='Make '+inputs['material_1']['formula']+' slab')
WF.append(Make_Slab_M1)

Make_Slab_M2 = Firework(FT_MakeSlabFromStructure(
                        bulk_name=inputs['material_2']['formula']+'bulk',
                        dict_name='material_2'),
                name='Make '+inputs['material_2']['formula']+' slab')
WF.append(Make_Slab_M2)

Relax_Slab_M1 = RelaxFW(name='Relax '+inputs['material_1']['formula']+' slab')
WF.append(Relax_Slab_M1)

Relax_Slab_M2 = RelaxFW(name='Relax '+inputs['material_2']['formula']+' slab')
WF.append(Relax_Slab_M2)

Make_Hetero_Structure = Firework(FT_MakeHeteroStructure(
    bottom_slab_name=inputs['material_1']['formula']+
                     inputs['material_1']['miller'],
    top_slab_name=inputs['material_2']['formula']+
                  inputs['material_2']['miller'],
    parameters_loc=['interface_params']), name='Make the interface')
WF.append(Make_Hetero_Structure)

Print_spec = Firework(FT_PrintSpec(), name='Print Spec at the end')
WF.append(Print_spec)

#Define dependencies:
Dependencies = {Initialize: [Check_Inputs],
                Check_Inputs: [Start_M1, Start_M2],
                Start_M1: [Converge_M1],
                Start_M2: [Converge_M2],
                Converge_M1: [Final_Params],
                Converge_M2: [Final_Params],
                Final_Params: [Relax_M1, Relax_M2],
                Relax_M1: [Make_Slab_M1],
                Relax_M2: [Make_Slab_M2],
                Make_Slab_M1: [Relax_Slab_M1],
                Make_Slab_M2: [Relax_Slab_M2],
                Relax_Slab_M1: [Make_Hetero_Structure],
                Relax_Slab_M2: [Make_Hetero_Structure],
                Make_Hetero_Structure: [Print_spec]}

wf = Workflow(WF, Dependencies, name='Dummy Heterogeneous Workflow')

# finally, instatiate the LaunchPad and add the workflow to it
lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(wf)

rapidfire(lpad)