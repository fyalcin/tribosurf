#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 11:28:20 2021

@author: yalcin
"""

from fireworks import LaunchPad, Workflow, Firework

from triboflow.firetasks.init_check import FTCheckInput, material_from_mp
from triboflow.firetasks.init_db import FT_PutMaterialInDB
from triboflow.firetasks.structure_manipulation import FT_StartPreRelax
from triboflow.workflows.subworkflows import surface_energy_swf
from triboflow.firetasks.dielectric import FT_StartDielectric
from triboflow.firetasks.convergence import FT_StartConvo
from triboflow.firetasks.init_check import unbundle_input

lpad = LaunchPad.auto_load()

db_file = 'auto'
high_level = True
override = True
fake = False

inputs = {'material': {'formula': 'ZnCu',
                       'miller': '111',
                       'mpid': 'mp-987',
                       'thick_min': 3,
                       'thick_max': 12,
                       'thick_incr': 1,
                       },
          'comp_params': {'functional': 'PBE',
                          'volume_tolerance': 0.001,
                          'BM_tolerance': 0.01,
                          'use_vdw': 'No',
                          'surfene_thr': 0.01,
                          'vacuum': 12},
          'sg_params': {'miller': [(0, 0, 1)],
                        'symmetrize': False,
                        'slab_thick': 12,
                        'vac_thick': 15,
                        'prim': True,
                        'lll_reduce': True,
                        'minimize_bv': True,
                        'max_index': 1,
                        'tol': 0.1,
                        'max_normal_search': 'max',
                        'reconstruct': True},
          'sg_filter': {'method': 'all'},
          'comp_params_user': {}
          }

material, comp_params, sg_params, sg_filter, comp_params_user = unbundle_input(inputs, keys=['material',
                                                                                             'comp_params',
                                                                                             'sg_params',
                                                                                             'sg_filter',
                                                                                             'comp_params_user'])

struct, mpid = material_from_mp(material)

functional = comp_params.get('functional', 'PBE')

ft_mat = FTCheckInput(input_dict=material,
                      read_key='material_params',
                      output_dict_name='material')

ft_computation = FTCheckInput(input_dict=comp_params,
                              read_key='comp_params',
                              output_dict_name='comp_params')

ft_mat_db = FT_PutMaterialInDB(mat='material', comp_params='comp_params')

InitFW = Firework([ft_mat, ft_computation, ft_mat_db])

PreRelaxation = Firework(FT_StartPreRelax(mp_id=mpid,
                                          functional=functional),
                         name='Start pre-relaxation for {}'
                         .format(material['formula']))

ConvergeEncut = Firework(FT_StartConvo(conv_type='encut',
                                       mp_id=mpid,
                                       functional=functional,
                                       ),
                         name='Start encut convergence for {}'
                         .format(material['formula']))

ConvergeKpoints = Firework(FT_StartConvo(conv_type='kpoints',
                                         mp_id=mpid,
                                         functional=functional,
                                         ),
                           name='Start kpoints convergence for {}'
                           .format(material['formula']))

CalcDielectric = Firework(FT_StartDielectric(mp_id=mpid,
                                             functional=functional,
                                             update_bulk=True,
                                             update_slabs=False),
                          name=f'Start dielectric SWF for {material["formula"]}')

WF1 = Workflow([InitFW, PreRelaxation, ConvergeEncut, ConvergeKpoints, CalcDielectric],
               {InitFW: [PreRelaxation],
                PreRelaxation: [ConvergeEncut],
                ConvergeEncut: [ConvergeKpoints],
                ConvergeKpoints: [CalcDielectric]})

WF2 = surface_energy_swf(mpid=mpid,
                         functional=functional,
                         sg_params=sg_params,
                         sg_filter=sg_filter,
                         db_file=db_file,
                         high_level=high_level,
                         comp_params_user=comp_params_user,
                         custom_id=None)

WF1.append_wf(WF2, WF1.leaf_fw_ids)

lpad.add_wf(WF1)
