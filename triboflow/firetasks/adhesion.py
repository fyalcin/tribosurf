#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 15:36:52 2021

@author: wolloch
"""
import numpy as np

from pymatgen.core.structure import Structure

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator

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
                    collection=nav.db.tasks,
                    fltr={'task_label': top_label})
        top_energy = top_calc['output']['energy']
        
        bot_calc = nav.find_data(
                    collection=nav.db.tasks,
                    fltr={'task_label': bot_label})
        bot_energy = bot_calc['output']['energy']
        
        inter_calc = nav.find_data(
                    collection=nav.db.tasks,
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