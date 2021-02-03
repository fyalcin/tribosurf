#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:09:03 2021

Firetasks to converge the slab thickness either calculating the surface 
energy or the lattice parameter for structures with increasing number of 
atomic layers.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction
from atomate.utils.utils import env_chk
from pymatgen.core.structure import Structure
from pymatgen.transformations.advanced_transformations import SlabTransformation

from triboflow.utils.database import Navigator, NavigatorMP
from triboflow.firetasks.slabs_wfs import SlabWFs, SlabThicknessError
from triboflow.tasks.io import read_json

currentdir = os.path.dirname(__file__)

# ============================================================================
# Firetasks
# ============================================================================

@explicit_serialize
class FT_StartSlabThicknessConvo(FiretaskBase):
    """
    Firetask description...
    
    """
    
    _fw_name = 'Slab Thickness convergence'

    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'high_level', 'relax_type', 'convo_kind'] #'bulk_struct_name', 'slab_out_name'
    
    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 
     
        # Read required parameters
        flag = self['mp_id']
        miller = self['miller']
        functional = self['functional']

        # Set optional parameters
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        high_level = self.get('high_level', None)
        convo_kind = self.get('convo_kind', 'surfene')

        defaults = read_json(currentdir + '/defaults.json')
        
        # Retrieve the bulk information from the high level DB
        nav = Navigator(db_file, high_level)
        slab = nav.find_data(functional+'.slab_data', {'mpid': flag, 
                                                       'miller': miller})

        # Run the Subworkflow to converge the thickness
        stop_convergence = slab.get('opt_thickness')
        if not stop_convergence:
            # Find out the structure
            bulk = nav.find_data(functional+'.bulk_data', {'mpid': flag})
            structure = Structure.from_dict(bulk.get('structure_fromMP'))
            comp_params = slab.get('comp_parameters', {})

            wf = select_thickness_convo(structure, flag, miller, comp_params, 
                                        functional, convo_kind)

            return FWAction(detours=wf, update_spec=fw_spec)

        # Continue the Workflow
        else:
            return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_OptimalThickness(FiretaskBase):

    """
    Firetask description...
    
    """
    
    _fw_name = 'Surface Energy calculation'
    required_params = ['flag', 'miller', 'functional', 'tag']
    optional_params = ['db_file', 'struct_out_name', 'file_output',
                       'output_dir', 'remote_copy', 'server', 'user', 'port']

    def run_task(self, fw_spec):
        pass

# ============================================================================
# Functions
# ============================================================================

def select_thickness_convo(structure, mp_id, miller, comp_params, 
                           functional='PBE', convo_kind='surfene'):
    """
    Select the kind of thickness convergence
    """
    
    if convo_kind == 'surfene':
        generate_wf = SlabWFs.converge_thickness_surfene
    elif convo_kind == 'alat':
        generate_wf = SlabWFs.converge_thickness_alat
    else:
        raise SlabThicknessError("Wrong input argument for convo_kind. "
                                 "Allowed options: 'surfene', 'alat'.")

    wf = generate_wf(structure = structure, flag = mp_id, miller = miller,
                      comp_parameters = comp_params, functional = functional, 
                      print_help = False)
    return wf
