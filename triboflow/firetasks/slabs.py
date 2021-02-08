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
    
    _fw_name = 'Converge the slab thickness'

    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'high_level', 'relax_type', 'convo_kind', 
                       'bulk_name', 'slab_name'] 
    
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
        relax_type = self.get('relax_type', 'slab_pos_relax')
        convo_kind = self.get('convo_kind', 'surfene')

        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        
        bulk_name = self.get('bulk_struct_name', 'structure_equiVol')
        slab_out_name = self.get('slab_out_name', 'relaxed_slab')

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
        
        bulk lattice = Lattice.cubic(3.508)
        bulk = Structure(lattice, ["Cu", "Cu", "Cu", "Cu"],
                    [[0,0,0], [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0]
                    ])
        miller = (1,1,0)
        min_thickness = 6
        min_vacuum = 0

        Potcar()

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

def read_slaboptparams():