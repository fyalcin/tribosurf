#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 15:36:52 2021

@author: wolloch
"""
import numpy as np
from uuid import uuid4

from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator
from triboflow.utils.structure_manipulation import (
    interface_name, slab_from_structure)
from triboflow.workflows.base import dynamic_relax_swf
from triboflow.utils.vasp_tools import get_custom_vasp_relax_settings


@explicit_serialize
class FT_RetrievMatchedSlabs(FiretaskBase):
    required_params = ['mp_id_1', 'mp_id_2', 'miller_1', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'top_out_name', 'bottom_out_name',
                       'high_level_db']
    
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import adhesion_energy_swf
        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        miller_1 = self.get('miller_1')
        miller_2 = self.get('miller_2')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        top_out_name = self.get('top_out_name', 'top_aligned_relaxed')
        bot_out_name = self.get('bottom_out_name', 'bottom_aligned_relaxed')
        hl_db = self.get('high_level_db', 'triboflow')
        
        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)
        
        input_list = fw_spec.get('relaxation_inputs')
        
        if input_list:
            nav = Navigator(db_file)
            for i in input_list:
                label = i[-1]
                miller = i[0].miller_index
                calc = nav.find_data(collection='tasks',
                                     fltr={'task_label': label})
                out_struct = calc['output']['structure']
                slab = slab_from_structure(miller,
                                           Structure.from_dict(out_struct))
                if label.startswith('top'):
                    out_name = top_out_name
                else:
                    out_name = bot_out_name
                nav_high = Navigator(db_file, high_level=hl_db)
                nav_high.update_data(collection=functional+'.interface_data',
                                     fltr={'name': name},
                                     new_values={'$set':
                                                 {out_name: slab.as_dict()}})
        
        return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_RelaxMatchedSlabs(FiretaskBase):
    required_params = ['mp_id_1', 'mp_id_2', 'miller_1', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'top_in_name', 'top_out_name',
                       'bottom_in_name', 'bottom_out_name', 'high_level_db']
    
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import adhesion_energy_swf
        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        miller_1 = self.get('miller_1')
        miller_2 = self.get('miller_2')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        top_in_name = self.get('top_in_name', 'top_aligned')
        top_out_name = self.get('top_out_name', 'top_aligned_relaxed')
        bot_in_name = self.get('bottom_in_name', 'bottom_aligned')
        bot_out_name = self.get('bottom_out_name', 'bottom_aligned_relaxed')
        hl_db = self.get('high_level_db', 'triboflow')
        
        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)

        nav = Navigator(db_file, high_level=hl_db)
        
        interface_dict = nav.find_data(collection=functional+'.interface_data',
                                       fltr={'name': name})
        
        relaxed_top_present = interface_dict.get(top_out_name)
        relaxed_bot_present = interface_dict.get(top_out_name)
        comp_params = interface_dict.get('comp_parameters', {})
        
        inputs = []
        if not relaxed_top_present:
            top_slab = Slab.from_dict(interface_dict[top_in_name])
            top_vis = get_custom_vasp_relax_settings(top_slab,
                                                     comp_params,
                                                     'slab_pos_relax')
            formula = top_slab.composition.reduced_formula
            miller = ''.join(str(s) for s in top_slab.miller_index)
            label = 'top_slab_'+formula+miller+'_'+str(uuid4())
            inputs.append([top_slab, top_vis, label])
        if not relaxed_bot_present:
            bot_slab = Slab.from_dict(interface_dict[bot_in_name])
            bot_vis = get_custom_vasp_relax_settings(bot_slab,
                                                     comp_params,
                                                     'slab_pos_relax')
            formula = bot_slab.composition.reduced_formula
            miller = ''.join(str(s) for s in bot_slab.miller_index)
            label = 'bot_slab_'+formula+miller+'_'+str(uuid4())
            inputs.append([bot_slab, bot_vis, label])
            
        if inputs:
            WF = dynamic_relax_swf(inputs, 'Relaxing the matched slabs')
            return FWAction(detours=WF,
                            update_spec={'relaxation_inputs': inputs})
        else:
            return FWAction(update_spec={'relaxation_inputs': inputs})
        
@explicit_serialize
class FT_StartAdhesionSWF(FiretaskBase):
    required_params = ['mp_id_1', 'mp_id_2', 'miller_1', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'adhesion_handle', 'high_level_db']
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import adhesion_energy_swf
        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        miller_1 = self.get('miller_1')
        miller_2 = self.get('miller_2')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        adhesion_handle = self.get('adhesion_handle', 'adhesion_energy@min')
        hl_db = self.get('high_level_db', 'triboflow')
            
        nav = Navigator(db_file, high_level=hl_db)
        
        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)
        
        interface_dict = nav.find_data(collection=functional+'.interface_data',
                                       fltr={'name': name})
        
        adhesion_was_calculated = interface_dict.get(adhesion_handle)
        comp_params = interface_dict.get('comp_parameters',{})
        
        if not adhesion_was_calculated:
            top_slab = Slab.from_dict(interface_dict['top_aligned_relaxed'])
            bottom_slab = Slab.from_dict(interface_dict['bottom_aligned_relaxed'])
            interface = Structure.from_dict(interface_dict['relaxed_structure@min']) 
        
            SWF = adhesion_energy_swf(top_slab,
                                  bottom_slab,
                                  interface,
                                  interface_name=None,
                                  functional='PBE',
                                  comp_parameters=comp_params)
            
            return FWAction(detours=SWF)
        else:
            return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_CalcAdhesion(FiretaskBase):
    
    required_params = ['interface_name', 'functional', 'top_label',
                       'bottom_label', 'interface_label']
    optional_params = ['db_file', 'out_name', 'high_level_db']

    def run_task(self, fw_spec):
        
        name = self.get('interface_name')
        functional = self.get('functional')
        top_label = self.get('top_label')
        bot_label = self.get('bottom_label')
        inter_label = self.get('interface_label')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        out_name = self.get('out_name', 'adhesion_energy@min')
        hl_db = self.get('high_level_db', 'triboflow')
        
        nav = Navigator(db_file=db_file)
        
        top_calc = nav.find_data(
                    collection='tasks',
                    fltr={'task_label': top_label})
        top_energy = top_calc['output']['energy']
        
        bot_calc = nav.find_data(
                    collection='tasks',
                    fltr={'task_label': bot_label})
        bot_energy = bot_calc['output']['energy']
        
        inter_calc = nav.find_data(
                    collection='tasks',
                    fltr={'task_label': inter_label})
        inter_energy = inter_calc['output']['energy']
        struct = Structure.from_dict(inter_calc['output']['structure'])
        
        area = np.linalg.norm(
            np.cross(struct.lattice.matrix[0],
                     struct.lattice.matrix[1])
            )
        
        E_abs = (top_energy + bot_energy) - inter_energy

        # Convert adhesion energz from eV/Angstrom^2 to J/m^2        
        E_Jm2 = 16.02176565 * E_abs / area
        
        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional+'.interface_data',
            fltr={'name': name},
            new_values={'$set': {out_name: E_Jm2}})
        
        return FWAction(update_spec=fw_spec)
