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

import os

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction
from atomate.utils.utils import env_chk
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.transformations.advanced_transformations import SlabTransformation

from triboflow.utils.database import Navigator, NavigatorMP
from triboflow.firetasks.slabs_wfs import SlabWFs
from triboflow.firetasks.errors import SlabThicknessError
from triboflow.tasks.io import read_json

currentdir = os.path.dirname(__file__)

# ============================================================================
# Firetasks
# ============================================================================

@explicit_serialize
class FT_SlabOptThick(FiretaskBase):
    """
    Start a subworkflow as a detour to calculate the optimal thickness generated
    from a provided bulk structure and with a certain orientation.
    The thickness can be converged either by evaluating how the surface energy
    or the lattice parameters changes with the number of atomic planes.

    bulk_name : str or None, optional
        Name of the bulk structure in the bulk database (material is
        identified by mp_id, but there might be different structures of the
        same material.) The default is None.
    slab_name : str or None, optional
        Name of the slab to be put in the DB, otherwise identified by mp_id 
        and miller index. The default is None.

    """
    
    _fw_name = 'Start a swf to converge slab thickness'

    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'low_level', 'high_level', 'relax_type', 
                       'convo_kind', 'bulk_name', 'slab_name'] 

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file = dfl, default_key="SlabOptThick")
        
        # Retrieve the bulk information from the high level DB
        nav = Navigator(p['db_file'], p['high_level'])
        slab = nav.find_data(p['functional']+'.slab_data', 
                             {'mpid': p['mp_id'], 'miller': p['miller']})

        # Start a subworkflow to converge the thickness if not already done
        stop_convergence = slab.get('opt_thickness', None)
        if not stop_convergence:
            # Retrieve the bulk structure
            bulk = nav.find_data(p['functional']+'.bulk_data', 
                                 {'mpid': p['mp_id']})
            structure = Structure.from_dict(bulk.get('structure_fromMP'))
            comp_params = slab.get('comp_parameters', {})

            wf = select_slabthick_conv(structure = structure, 
                                       mp_id = p['mp_id'], 
                                       miller = p['miller'],
                                       functional = p['functional'], 
                                       comp_params = comp_params,
                                       db_file = p['db_file'],
                                       low_level = p['low_level'],
                                       high_level = p['high_level'],
                                       bulk_name = p['bulk_name'],
                                       slab_name = p['slab_name'])

            return FWAction(detours=wf, update_spec=fw_spec)

        # Continue the Workflow
        else:
            return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_StartThickConvo(FiretaskBase):
    """
    Firetask description...
    
    """
    
    _fw_name = 'Surface Energy calculation'
    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'low_level', 'high_level', 'relax_type', 
                       'convo_kind', 'bulk_name', 'slab_name', ] 

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """
        
        min_thickness = 6
        min_vacuum = 0

        slabgen = SlabGenerator(initial_structure = bulk,
                                miller_index = miller,
                                center_slab = True,
                                primitive = False,
                                lll_reduce = True,
                                in_unit_planes = True,
                                #max_normal_search=max([abs(l) for l in miller]),
                                min_slab_size = min_thickness,
                                min_vacuum_size = min_vacuum)

        slabgen = SlabGenerator(initial_structure = bulk,
                                miller_index = miller,
                                center_slab=True,
                                primitive=False,
                                lll_reduce=True,
                                in_unit_planes=True,
                                #max_normal_search=max([abs(l) for l in miller]),
                                min_slab_size=min_thickness,
                                min_vacuum_size=min_vacuum)
        bulk = slabgen.get_slab()

        #bulk.to(fmt='poscar', filename='POSCAR')
        slab.to(fmt='poscar', filename='POSCAR')

@explicit_serialize
class FT_EndThickConvo(FiretaskBase):
    """
    Firetask description...
    
    """

    def run_task(elf, fw_spec):
        """ Run the Firetask.
        """ 
    pass

@explicit_serialize
class FT_GenerateSlabs(FiretaskBase):
    """
    Firetask description...
    
    """

    def run_task(elf, fw_spec):
        """ Run the Firetask.
        """ 
    pass


# ============================================================================
# Functions
# ============================================================================

def select_slabthick_conv(structure, mp_id, miller, comp_params, 
                          functional='PBE', convo_kind='surfene'):
    """
    [summary]

    Parameters
    ----------
    structure : [type]
        [description]
    mp_id : [type]
        [description]
    miller : [type]
        [description]
    comp_params : [type]
        [description]
    functional : str, optional
        [description], by default 'PBE'
    convo_kind : str, optional
        [description], by default 'surfene'

    Returns
    -------
    [type]
        [description]

    Raises
    ------
    SlabThicknessError
        [description]
    """

    if convo_kind == 'surfene':
        generate_wf = SlabWFs.conv_slabthick_surfene
    elif convo_kind == 'alat':
        generate_wf = SlabWFs.conv_slabthick_alat
    else:
        raise SlabThicknessError("Wrong input argument for convo_kind. "
                                 "Allowed options: 'surfene', 'alat'.")

    wf = generate_wf(structure = structure, flag = mp_id, miller = miller,
                      comp_parameters = comp_params, functional = functional, 
                      print_help = False)
    return wf

def read_runtask_params(obj, fw_spec, required_params, optional_params,
                        default_file, default_key):
    """
    [summary]

    Parameters
    ----------
    obj : [type]
        [description]
    fw_spec : [type]
        [description]
    required_params : [type]
        [description]
    optional_params : [type]
        [description]
    default_file : [type]
        [description]
    default_key : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    defaults = read_json(default_file)
    params = {}

    # Read required and optional parameters
    for key in required_params:
        params[key] = obj[key]
    for key in optional_params:
        params[key] = obj.get(key, defaults[default_key][key])

    if not 'db_file' in params.keys():
        params['db_file'] = env_chk('>>db_file<<', fw_spec)

    return params


# # Read required parameters
# flag = self['mp_id']
# miller = self['miller']
# functional = self['functional']

# # Set optional parameters
# db_file = self.get('db_file')
# if not db_file:
#     db_file = env_chk('>>db_file<<', fw_spec)
# high_level = self.get('high_level', None)
# relax_type = self.get('relax_type', 'slab_pos_relax')
# convo_kind = self.get('convo_kind', 'surfene')

# slab_name = self.get('slab_struct_name', 'unrelaxed_slab')

# bulk_name = self.get('bulk_struct_name', 'structure_equiVol')
# slab_out_name = self.get('slab_out_name', 'relaxed_slab')