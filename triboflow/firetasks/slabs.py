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
import monty
from uuid import uuid4

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction
from fireworks import Workflow, Firework
from atomate.utils.utils import env_chk
from atomate.vasp.powerups import add_modify_incar
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.transformations.advanced_transformations import SlabTransformation

from triboflow.utils.database import StructureNavigator, NavigatorMP
from triboflow.firetasks.slabs_wfs import SlabWFs
from triboflow.firetasks.surfene_wfs import SurfEneWfs
from triboflow.utils.errors import SlabOptThickError, GenerateSlabsError, RelaxStructureError
from triboflow.tasks.io import read_json
from triboflow.utils.vasp_tools import GetCustomVaspRelaxSettings
from atomate.vasp.fireworks.core import OptimizeFW, ScanOptimizeFW

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
                       'convo_kind', 'thick_start', 'thick_incr', 'nsteps',
                       'vacuum', 'bulk_name', 'slab_name',]

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file = dfl, default_key="SlabOptThick")
        
        # Retrieve the bulk information from the high level DB
        nav_struct = StructureNavigator(p['db_file'], p['high_level'])
        slab = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                           functional=p['functional'], 
                                           miller=p['miller'])

        # Start a subworkflow to converge the thickness if not already done
        stop_convergence = slab.get('opt_thickness', None)
        if not stop_convergence:
            # Retrieve the bulk structure
            bulk = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                               functional=p['functional'])
            if p['bulk_name'] is not None:
                bulk = bulk[p['bulk_name']]

            structure = Structure.from_dict(bulk.get('structure_fromMP'))
            comp_params = slab.get('comp_parameters', {})

            wf = self.select_slabthick_conv(structure=structure, 
                                            mp_id=p['mp_id'], 
                                            miller=p['miller'],
                                            functional=p['functional'], 
                                            comp_params=comp_params,
                                            db_file=p['db_file'],
                                            low_level=p['low_level'],
                                            high_level=p['high_level'],
                                            thick_start=p['thick_start'],
                                            thick_incr=p['thick_incr'],
                                            nsteps=p['nsteps'],
                                            vacuum=p['vacuum'],
                                            slab_name=p['slab_name'],
                                            convo_kind=p['convo_kind'],
                                            relax_type=p['relax_type'])

            return FWAction(detours=wf, update_spec=fw_spec)

        # Continue the Workflow
        else:
            return FWAction(update_spec=fw_spec)

    def select_slabthick_conv(self, structure, mp_id, miller, comp_params, 
                              db_file=None, low_level=None, high_level='triboflow', 
                              thick_start=4, thick_incr=1, nsteps=6, vacuum=10, 
                              slab_name=None, functional='PBE', 
                              relax_type="slab_pos_relax", convo_kind='surfene'):
        """
        Select the desired subworkflow from the SlabWFs class, to converge the 
        slab thickness either via surface energy or lattice parameter.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            Structure of the bulk to be processed.

        mp_id : str
            [description]

        miller : str
            [description]

        comp_params : dict
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
            When a wrong input argument is passed to convo_kind
        """

        if convo_kind == 'surfene':
            generate_wf = SlabWFs.conv_slabthick_surfene
        elif convo_kind == 'alat':
            generate_wf = SlabWFs.conv_slabthick_alat
        else:
            raise SlabOptThickError("Wrong input argument for convo_kind. "
                                    "Allowed options: 'surfene', 'alat'.")

        wf = generate_wf(structure = structure, mp_id = mp_id, miller = miller,
                        functional = functional, comp_params = comp_params,
                        db_file=db_file, low_level=low_level, high_level=high_level,
                        thick_start=thick_start, thick_incr=thick_incr, nsteps=nsteps,
                        vacuum=vacuum, slab_name=slab_name, relax_type=relax_type)

        return wf

@explicit_serialize
class FT_StartThickConvo(FiretaskBase):
    """
    Firetask description...
    
    """
    
    _fw_name = 'Start the slab thickness convergence'
    required_params = ['structure', 'mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'collection', 'comp_params', 'convo_kind',
                       'relax_type', 'thick_start', 'thick_incr', 'nsteps',
                       'vacuum', 'ext_index', 'slab_name', 'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file = dfl, default_key="StartThickConvo")

        wf = SurfEneWfs.surface_energy_wf(structure=p['structure'], 
                                          mp_id=p['mp_id'], 
                                          miller=p['miller'], 
                                          functional=p['functional'], 
                                          db_file=p['db_file'], 
                                          collection=p['collection'],
                                          thick_start=p['thick_start'], 
                                          thick_incr=p['thick_incr'], 
                                          nsteps=p['nsteps'], 
                                          vacuum=p['vacuum'], 
                                          ext_index=p['ext_index'], 
                                          slab_name=p['slab_name'])
        
        return FWAction(detours=wf, update_spec=fw_spec)

@explicit_serialize
class FT_EndThickConvo(FiretaskBase):
    """
    Firetask description...
    
    """

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 
    pass

@explicit_serialize
class FT_GenerateSlabs(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['structure', 'mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'collection', 'thickness', 'vacuum', 
                       'ext_index', 'symmetrize', 'slab_name']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file = dfl, default_key="GenerateSlabs")

        GenerateSlabsError.check_slabname_thickness(p['thickness'], 
                                                    p['slab_name'])

        # Generate the slabs for each structure passed as input
        slabs = []
        thickness = list(p['thickness'])
        slab_name = list(p['slab_name'])

        for thk in thickness:
            slabgen = SlabGenerator(initial_structure = p['structure'],
                                    miller_index = p['miller'],
                                    center_slab=True,
                                    primitive=False,
                                    lll_reduce=True,
                                    in_unit_planes=True,  # Fundamental
                                    min_slab_size=thk,
                                    min_vacuum_size=p['vacuum'])
            s = slabgen.get_slabs(bonds=None, ftol=0.1, tol=0.1, 
                                  max_broken_bonds=0, symmetrize=p['symmetrize'],
                                  repair=False)
            # If exit index is True (default), return the first element.
            # TODO: Generalize this approach by evaluating the different slabs
            # that are returned and select the optimal one.
            if p['ext_index']:
                s = s[0]

            slabs.append(s)

        # Add the slabs to the "collection" Database
        nav_struct = StructureNavigator(p['db_file'], p['collection'])
        
        # Store unrelaxed data in the Database
        for slab, name in zip(slabs, slab_name):
            slab_dict = monty.json.jsanitize(s.as_dict(), allow_bson=True)
            nav_struct.add_slab_to_db(structure=slab, 
                                      mp_id=p['mp_id'],
                                      functional=p['functional'],
                                      miller=p['miller'],
                                      struct_name=name)
        
        # coll.update_one({'mpid': flag, 'miller': miller},
        #                 {'$set': {'unrelaxed_slab': slab_dict}},
        #                 upsert=True)
        
        return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_RelaxStructure(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['mp_id', 'functional', 'struct_kind',]
    optional_params = ['comp_params', 'miller', 'name', 'db_file', 'collection',
                       'relax_type']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file = dfl, default_key="RelaxStructure")

        # Get the structure from the Database with a calculation tag
        structure, tag = self.extract_data_from_db()
        if p['name'] is not None:
            structure = structure['name']

        # Set the simulation if the calculation is not already at convergence
        if 'relaxed' or 'relaxed_slab' not in structure:
            fw = self.set_calculation(structure, tag, dfl)

            # Set the Relaxation Workflow for making a detour.
            # Use add_modify_incar powerup to add KPAR and NCORE settings based on 
            # env_chk in my_fworker.yaml.
            wf_name = tag.split('_')[0] + p['relax_type']
            opt_wf = add_modify_incar(Workflow([fw], name=wf_name))

            return FWAction(detours=opt_wf)
    
    def extract_data_from_db(self):

        RelaxStructureError.check_struct_kind(self.p['struct_kind'])

        # Call the navigator for retrieving structure
        nav_struct = StructureNavigator(db_file=self.p['db_file'],
                                        high_level=self.p['high_level'])

        # Extract data from database
        if self.p['struct_kind'] == 'bulk':
            structure = nav_struct.get_bulk_from_db(self.p['mp_id'], 
                                                    self.p['functional'],
                                                    warning=True)
        elif self.p['struct_kind'] == 'slab':
            structure = nav_struct.get_slab_from_db(self.p['mp_id'], 
                                                    self.p['functional'], 
                                                    self.p['miller'], 
                                                    warning=True)
        elif self.p['struct_kind'] == 'interface':
            structure = nav_struct.get_interface_from_db(self.p['name'],
                                                         self.p['functional'],
                                                         warning=True)
        
        RelaxStructureError.is_data(structure=structure, mp_id=self.p['mp_id'], 
                                    functional=self.p['functional'], 
                                    struct_kind=self.p['struct_kind'])
        
        # Get a tag for the calculation
        formula = structure.composition.reduced_formula
        if self.p['miller'] is not None:
            tag = formula + str(p['miller']) + '_' + str(uuid4())
        else:
            tag = formula + '_' + str(uuid4())

        return structure, tag

    def set_calculation(self, structure, tag, dfl):

        # Check the computational parameters
        comp_params = self.p['comp_params']
        if not comp_params:
            defaults = read_json(dfl)
            comp_params = defaults['comp_params']

        # Set options for vasp
        vis = GetCustomVaspRelaxSettings(structure, comp_params,
                                         self.p['relax_type'])
        
        # Create the Firework to run the simulation
        if self.p['functional'] == 'SCAN':
            fw = ScanOptimizeFW(structure=structure, name=tag, vasp_input_set=vis)
        else:
            fw = OptimizeFW(structure, name=tag, vasp_input_set=vis,
                            half_kpts_first_relax=True)

        return fw

# ============================================================================
# Functions
# ============================================================================

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
