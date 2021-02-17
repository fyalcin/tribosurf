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
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

import os
import monty
from uuid import uuid4
from pprint import pprint, pformat

import numpy as np
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction
from fireworks import Firework, Workflow, FileWriteTask
from atomate.utils.utils import env_chk
from atomate.vasp.powerups import add_modify_incar
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, Slab
from atomate.vasp.fireworks.core import OptimizeFW, ScanOptimizeFW

from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.errors import (
    SlabOptThickError, 
    GenerateSlabsError, 
    RelaxStructureError,
    MoveStructInDBError,
    ReadParamsError,
    WriteParamsError
)
from triboflow.tasks.io import read_json
from triboflow.utils.vasp_tools import GetCustomVaspRelaxSettings
from triboflow.firetasks.surfene_wfs import get_miller_str
from triboflow.firetasks.slabs_wfs import read_default_params
from triboflow.utils.file_manipulation import copy_output_files

currentdir = os.path.dirname(__file__)

# ============================================================================
# Firetasks
# ============================================================================

@explicit_serialize
class FT_SlabOptThick(FiretaskBase):
    """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

    Start a subworkflow as a detour to calculate the optimal thickness generated
    from a provided bulk structure and with a certain orientation.
    The thickness can be converged either by evaluating how the surface energy
    or the lattice parameters changes with the number of atomic planes.
    At the moment only the surface energy convergence is implemented.

    Parameters
    ----------
    mp_id : str
        MP-id of the structure from the MP database.
    
    miller : list of int, or str
        Miller indexes (h, k, l) to select the slab orientation.
        
    functional : str
        Functional for the pseudopotential to be adopted.

    db_file : str or None
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.
        
    low_level : str or None, optional
        Name of the table of the "low level" database, saved in db_file. The 
        intermediate calculations and raw data during will be saved here. If 
        nothing is passed the Firework database is used. The default is None.

    high_level : str, optional
        Name of the table of the "high level" database, saved in db_file.
        The slab optimal thickness will be saved here. The slab energy and
        surface energy will be saved too if convo_kind is 'surfene'.
        The default is 'triboflow'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.

    convo_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene', 
        'alat'. The latter is not implemented yet. The default is 'surfene'.
    
    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point.
        The default is 4.
        
    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If the convergence
        is not reached this value is assumed to be the optimal one.
        The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step.
        The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.
    
    ext_index : int, optional
        Use the ext_index element from SlabGenerator.get_slabs as a slab.
        The default is 0.

    bulk_name : str or list None, optional
        Name of the custom bulk dictionary, to be retrieved from the high level
        database and to be used to build the slabs. Bulks are identified by 
        mp_id and functional but there might be different structures of the
        same material. The default is "structure_fromMP".

    slab_name : str or None, optional
        Where to search for the information about the optimal thickness in the
        slab dictionary in the high level database. The default is None.

    """
    
    _fw_name = 'Start a subworkflow to converge slab thickness'

    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'low_level', 'high_level', 'convo_kind',
                       'relax_type', 'thick_min', 'thick_max', 'thick_incr',
                       'vacuum', 'in_unit_planes', 'ext_index', 'bulk_name', 
                       'slab_name']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, self.optional_params,
                                default_file=dfl, default_key="SlabOptThick")
        
        # Retrieve the bulk information from the high level DB
        nav_struct = StructureNavigator(p['db_file'], p['high_level'])
        slab = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                           functional=p['functional'], 
                                           miller=p['miller'])
        
        # If data is saved elsewhere from standard position
        if p['slab_name'] is not None:
            slab = get_one_info_from_struct_dict(slab, p['slab_name'])

        # Start a subworkflow to converge the thickness if not already done
        stop_convergence = slab.get('opt_thickness', None)

        if not stop_convergence:
            
            # Retrieve the desired bulk structure
            bulk = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                               functional=p['functional'])
            bulk = get_one_info_from_struct_dict(bulk, p['bulk_name'])

            structure = Structure.from_dict(bulk.get('structure_fromMP'))
            comp_params = slab.get('comp_params', {})

            wf = self.select_slabthick_conv(structure=structure, 
                                            comp_params=comp_params,
                                            p=p)

            return FWAction(detours=wf, update_spec=fw_spec)

        # Continue the Workflow
        else:
            return FWAction(update_spec=fw_spec)

    def select_slabthick_conv(self, structure, comp_params, p):
        """
        Select the desired subworkflow from the SlabWFs class, to converge the 
        slab thickness either by evaluating the surface energy or the lattice 
        parameter.
        
        """

        from triboflow.firetasks.slabs_wfs import SlabWF

        if p['convo_kind'] == 'surfene':
            generate_wf = SlabWF.conv_slabthick_surfene
        elif p['convo_kind'] == 'alat':
            generate_wf = SlabWF.conv_slabthick_alat
        else:
            raise SlabOptThickError("Wrong input argument for convo_kind. "
                                    "Allowed options: 'surfene', 'alat'")

        wf = generate_wf(structure=structure, mp_id=p['mp_id'], 
                         miller=p['miller'], functional=p['functional'], 
                         comp_params=p['comp_params'], db_file=p['db_file'],
                         low_level=p['low_level'], high_level=p['high_level'],
                         relax_type=p['relax_type'], thick_min=p['thick_min'], 
                         thick_max=p['thick_max'], thick_incr=p['thick_incr'], 
                         vacuum=p['vacuum'], in_unit_planes=p['in_unit_planes'],
                         ext_index=p['ext_index'])

        return wf

@explicit_serialize
class FT_StartThickConvo(FiretaskBase):
    """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

    It starts a subworkflow as a detour to converge the thickness of a slab.
    The thickness can be converged either by evaluating how the surface energy
    or the lattice parameters changes with the number of atomic planes.
    It is the first element of the SlabWF.conv_slabthick_surfene workflow.

    (At the moment only the surface energy convergence is implemented)
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Bulk structure to build the slabs.

    mp_id : str
        MP-id of the structure from the MP database.
    
    miller : list of int, or str
        Miller indexes (h, k, l) to select the slab orientation.
        
    functional : str
        Functional for the pseudopotential to be adopted.

    db_file : str or None
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.
        
    database : str or None, optional
        Name of the table where the data will be stored, within db_file. If 
        nothing is passed the Firework database is used. The default is None.

    convo_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene', 
        'alat'. The latter is not implemented yet. The default is 'surfene'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.

    comp_params : dict, optional
        Computational parameters for the VASP simulations. If not set, default 
        parameters will be used instead. The default is {}.
    
    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point.
        The default is 4.
        
    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If the convergence
        is not reached this value is assumed to be the optimal one.
        The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step.
        The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.

    ext_index : int, optional
        Use the ext_index element from SlabGenerator.get_slabs as a slab.
        The default is 0.

    cluster_params : dict, optional
        Dictionary containing cluster-related options to run efficiently the
        VASP simulations on a cluster. The default is {}.

    """
    
    _fw_name = 'Start the slab thickness convergence'
    required_params = ['structure', 'mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'low_level', 'high_level', 'convo_kind', 
                       'relax_type', 'comp_params', 'thick_min', 'thick_max', 
                       'thick_incr', 'vacuum', 'in_unit_planes', 'ext_index', 
                       'parallelization', 'recursion', 'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """
        from triboflow.firetasks.surfene_wfs import SurfEneWF

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, 
                                self.optional_params, default_file=dfl, 
                                default_key="StartThickConvo")
        
        # Select the convergence of interest
        wf = self.select_conv(p)

        return FWAction(detours=wf, update_spec=fw_spec)

        def select_conv(self, p):

            if p['convo_kind'] == 'surfene':
                wf = SurfEneWF.surface_energy(structure=p['structure'],
                                            mp_id=p['mp_id'], 
                                            miller=p['miller'], 
                                            functional=p['functional'],
                                            db_file=p['db_file'], 
                                            low_level=p['low_level'],
                                            high_level=p['high_level'],
                                            relax_type=p['relax_type'],
                                            comp_params=p['comp_params'],
                                            thick_min=p['thick_min'], 
                                            thick_max=p['thick_max'],
                                            thick_incr=p['thick_incr'],
                                            vacuum=p['vacuum'],
                                            in_unit_planes=p['slab_name'],
                                            ext_index=p['ext_index'], 
                                            parallelization=p['parallelization'],
                                            recursion=p['recursion'],
                                            cluster_params=p['cluster_params'])

            else:
                raise SystemExit('Lattice parameter convergence not yet implemented')


@explicit_serialize
class FT_EndThickConvo(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['structure', 'mp_id', 'miller']
    optional_params = ['db_file', 'low_level', 'high_level', 'functional', 
                       'convo_kind', 'relax_type', 'comp_params', 'thick_min', 
                       'thick_max', 'thick_incr', 'vacuum', 'in_unit_planes', 
                       'ext_index', 'parallelization', 'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, 
                                self.optional_params, default_file=dfl, 
                                default_key="EndThickConvo")
        
        # Retrieve the surface energy
        surface_energies = self.get_data(p)

        # Analyze energies, if convergence is reached data is stored in DB
        stop_convergence = self.analyze_data(surface_energies, p)

        # Decide whether to rerun recursively the surface energy workflow
        if not stop_convergence:
            wf = None
            return FWAction(detours=wf, update_spec=fw_spec)

        else:
            return FWAction(update_spec=fw_spec)
        
        def get_data(self, p):

            if p['convo_kind'] == 'surfene':
                data = self._get_surfene(p)

            else:
                raise SlabOptThickError("Wrong argument: 'convo_kind'. Allowed "
                                        "values: 'surfene', 'alat'")
            
            return data
        
        def analyze_data(self, energies, p):

            if p['convo_kind'] == 'surfene':
                data = self._analyze_surfene(energies, p)
                
            else:
                raise SlabOptThickError("Wrong argument: 'convo_kind'. Allowed "
                                        "values: 'surfene', 'alat'")
            
            return stop_convergence
        
        # def recursion_wf(self, p):
        #     wf = generate_wf(structure=structure, mp_id=p['mp_id'], 
        #             miller=p['miller'], functional=p['functional'], 
        #             comp_params=p['comp_params'], db_file=p['db_file'],
        #             low_level=p['low_level'], high_level=p['high_level'],
        #             relax_type=p['relax_type'], thick_min=p['thick_min'], 
        #             thick_max=p['thick_max'], thick_incr=p['thick_incr'], 
        #             vacuum=p['vacuum'], in_unit_planes=p['in_unit_planes'],
        #             ext_index=p['ext_index']

        def _get_surfene(self, p):
            
            # Create the list of keys to be used when reading the dictionary
            thickness = np.arange(p['thick_min'], p['thick_max'], p['thick_incr'])
            name = [['thickness', 'data_' + str(thk), 'calc_output', 'surface_energy'] 
                    for thk in thickness]

            # Call the navigator for retrieving the dictionary out of the DB
            nav = Navigator(db_file=p['db_file'], high_level=p['high_level'])
            filter = {'mpid': p['mp_id'], 'miller': p['miller']}
            data_dict = nav.find_data(p['functional'] + '.slab_data', filter)

            # Get the surface energies
            energies = get_multiple_info_from_struct_dict(data_dict, p['name'])

            return energies
        
        def analyze_surfene(self, energies, p):
            pass


@explicit_serialize
class FT_GenerateSlabs(FiretaskBase):
    """
     Generate a slab or a list of slabs out of a given structure. Parameters 
     hat are taken into account to generate the possible different slabs are: 
     miller, thickness, vacuum, slab_name. The slabs are generated with 
     SlabGenerator and stored in the database.
    
    """

    required_params = ['structure', 'mp_id', 'miller', 'collection']
    optional_params = ['db_file', 'database', 'thickness', 'thick_max', 'vacuum',
                       'symmetrize', 'ext_index', 'in_unit_planes', 'slab_name']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, self.optional_params,
                                default_file=dfl, default_key="GenerateSlabs")

        # Generate the slabs for each structure passed as input
        miller, thickness, vacuum, slab_name = self.set_lists_to_loop(p)

        # Generate the slabs, taking into account miller, thickness, vacuum
        slabs = generate_slabs(structure=p['structure'], 
                               miller=miller, 
                               thickness=thickness, 
                               vacuum=vacuum, 
                               thick_bulk=p['thick_max'],
                               ext_index=p['ext_index'],
                               in_unit_planes=p['in_unit_planes'])
        
        # Store the slab structure in collection.name in database, within db_file
        self.structure_in_db(slabs, miller, slab_name, p)
        
        return FWAction(update_spec=fw_spec)
    
    def set_lists_to_loop(self, p):

        GenerateSlabsError.check_inputs(p['miller'], p['thickness'], 
                                        p['vacuum'], p['slab_name'])

        miller = p['miller']
        thickness = list(p['thickness'])
        vacuum = p['vacuum']
        slab_name = list(p['slab_name'])
        
        if not all([isinstance(x, list) for x in miller]):
            miller = [miller] * len(thickness)
        
        if not isinstance(vacuum, list):
            vacuum = [vacuum] * len(thickness)

        return miller, thickness, vacuum, slab_name
    
    def structure_in_db(self, slabs, miller, slab_name, p):

        nav = Navigator(p['db_file'], p['database'])
        
        # Store unrelaxed data in the Database
        for s, hkl, name in zip(slabs, miller, slab_name):
            # Clean the data and create a dictionary with the given path
            slab_dict = monty.json.jsanitize(s.as_dict(), allow_bson=True)
            update_data = write_one_dict_for_db(slab_dict, name)

            nav.update_data(collection=p['collection'], 
                            filter={'mpid': p['mp_id'], 'miller': hkl},
                            new_values={'$set': update_data},
                            upsert=True)


@explicit_serialize
class FT_RelaxStructure(FiretaskBase):
    """
    Retrieve a structure based on mp_id, and name out of a collection found in
    db_file and database. The slab is relaxed following the 'relax_type' procedure
    using an Atomate workflow. The result is stored in the same database with 
    a tag that can be provided by the user.
    
    """

    required_params = ['mp_id', 'functional', 'collection', 'name', 'tag']
    optional_params = ['db_file', 'database', 'relax_type', 'comp_params', 
                       'miller', 'check_key']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self,
                                fw_spec, 
                                self.required_params, 
                                self.optional_params,
                                default_file=dfl, 
                                default_key="RelaxStructure")

        # Retrieve the structure and check if it has been already calculated
        structure, is_done = getcheck_struct(p, pymatgen_obj=False)

        # Run the simulation if the calculation is not already done
        if not is_done:
            wf = self.run_relax_detour(structure, p, dfl)
            return FWAction(detours=wf, update_spec=fw_spec)

        # Continue the Workflow
        else:
            return FWAction(update_spec=fw_spec)

    def getcheck_struct(self, p, pymatgen_obj=False):
        
        # Check if collection does exist
        RelaxStructureError.check_collection(p['collection'])

        # Retrieve the structure from the Database
        structure = retrieve_from_db(db_file=p['db_file'], 
                                     database=p['database'], 
                                     collection=p['collection'], 
                                     mp_id=p['mp_id'],
                                     miller=p['miller'],
                                     name=p['name'],
                                     pymatgen_obj=pymatgen_obj)
        RelaxStructureError.is_data(structure, p['mp_id'], p['functional'])

        # Check if the calculation is already done, searching for given keys
        if p['check_key'] is not None:
            is_done = True if p['check_key'] in structure.keys() else False

        return structure, is_done      

    def run_relax_detour(self, structure, p, dfl):

        # Check tag and computational parameters
        tag = p['tag']
        comp_params = p['comp_params']
        if not bool(comp_params):
            comp_params = read_default_params(dfl, 'comp_params', comp_params)
        
        # Set options for vasp
        vis = GetCustomVaspRelaxSettings(structure, comp_params, p['relax_type'])
        
        # Create the Firework to run perform the simulation
        if p['functional'] == 'SCAN':
            fw = ScanOptimizeFW(structure=structure, name=tag, vasp_input_set=vis)
        else:
            fw = OptimizeFW(structure, name=tag, vasp_input_set=vis,
                            half_kpts_first_relax=True)

        # Define the workflow name
        wf_name = p['mp_id'] + '_' + p['relax_type']

        # Use add_modify_incar to add KPAR and NCORE settings based on env_chk
        wf = add_modify_incar(Workflow([fw], name=wf_name))

        return wf


@explicit_serialize
class FT_MoveStructInDB(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['mp_id', 'collection_from', 'collection_to', 'tag']
    optional_params = ['db_file', 'database_from', 'database_to', 'miller',
                       'name_check', 'name', 'name_tag', 'struct_kind', 
                       'override', 'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, 
                                fw_spec,
                                self.required_params,
                                self.optional_params,
                                default_file=dfl,
                                default_key="MoveStructInDB")

        # Check if a structure is already present in name and check_key
        is_done = self.check_struct(p)
        
        if not is_done or (is_done and p['override']):
            
            # Extract the information and store in destination db
            vasp_calc, info = self.get_results_from_tag(p)
            self.store_results(info, p)

            # Manage stdout, save a local poscar with results
            wf = self.user_output(vasp_calc, p)
            if wf is not None:
                return FWAction(detours=wf, update_spec=fw_spec)
        
        else:  
            return FWAction(update_spec=fw_spec)
    
    def check_struct(self, p):
        
        # Check if collection does exist
        MoveStructInDBError.check_collection(p['collection'])

        if p['name_check'] is None:
            is_done = False
        else:
            # Retrieve the structure from the Database
            structure = retrieve_from_db(db_file=p['db_file'], 
                                         database=p['database'], 
                                         collection=p['collection'], 
                                         mp_id=p['mp_id'],
                                         miller=p['miller'],
                                         name=p['name_check'],
                                         pymatgen_obj=False)
        
            # Check if the calculation is already done
            is_done = False if (structure is None or not bool(structure)) else False

        return is_done
    
    def get_results_from_tag(self, p):
        
        # Retrieve the vasp_calc_output and the info
        vasp_calc, info = retrieve_from_tag(db_file=p['db_file'],
                                            collection=p['collection_from'],
                                            tag=p['tag'],
                                            name=p['name_tag'],
                                            database=p['database_from'])
        return vasp_calc, info
    
    def store_results(self, info, p):

        # Dictionaries are stored in name[i]/name_tag[i]
        name = []
        for n, n_t in zip(list(p['name']), list(p['name_tag'])):
            name.append(list(n).append(list(n_t)[-1]))
        
        # Prepare the list of dictionaries to be stored in the database
        info_dict = write_multiple_dict_for_db(info, name)
    
        # Prepare the database and options where to store data
        nav = Navigator(db_file=p['db_file'], high_level=p['database_to'])
        filter = {'mpid': p['mp_id']}
        if p['miller'] is not None:
            filter.update({'miller': p['miller']})
    
        # Effectively store the data
        for d in info_dict:
            nav.update_data(p['collection_to'], filter, {'$set': d})
    
    def user_output(self, vasp_calc, p, dfl):

        cluster_params = p['cluster_params']
        
        # Set missing values in cluster parameters
        cluster_params = read_default_params(default_file=dfl, 
                                             default_key="cluster_params", 
                                             cluster_params=cluster_params)

        func = _select_func_from_dict(p['struct_kind'])
        structure = func.from_dict(vasp_calc['output']['structure'])

        # Output to screen
        print('')
        print('Relaxed output structure as pymatgen.surface.Slab dictionary:')
        pprint(structure.as_dict())
        print('')
        
        # handle file output:
        if p['file_output']:
            
            # Define POSCAR and dictionary ouput names
            prefix = p['mp_id'] + p['struct_kind']
            if p['miller'] is not None:
                prefix = prefix + p['miller']
            poscar_name = prefix + '_POSCAR.vasp'
            structure_name = prefix + '_dict.txt'

            # Define the subworkflow to write and copy the structure
            poscar_str = Poscar(structure).get_string()
            write_ft = FileWriteTask(files_to_write=
                                     [{'filename': poscar_name,
                                       'contents': poscar_str},
                                      {'filename': structure_name,
                                       'contents': pformat(structure.as_dict())}])
            copy_ft = copy_output_files(file_list=[poscar_name, structure_name],
                                        output_dir=p['output_dir'],
                                        remote_copy=p['remote_copy'],
                                        server=p['server'],
                                        user=['server'],
                                        port=['port'])
            fw = Firework([write_ft, copy_ft],
                          name='Copy the results of structure result')
            wf = Workflow.from_Firework(fw, name='Copy structure to file')

        else:
            wf = None
        
        return wf


# ============================================================================
# Functions
# ============================================================================

def read_runtask_params(obj, fw_spec, required_params, optional_params,
                        default_file, default_key):
    """
    General function to read the required and optional parameters of a 
    firetask instance. If 

    """

    defaults = read_json(default_file)
    params = {}

    # Read required and optional parameters
    for key in required_params:
        params[key] = obj[key]
    for key in optional_params:
        params[key] = obj.get(key, defaults[default_key][key])

    # Clean db_file
    if 'db_file' in params.keys():
        if params['db_file'] is None:
            params['db_file'] = env_chk('>>db_file<<', fw_spec)

    # Clean miller index
    if 'miller' in params.keys():
        if isinstance(params['miller'], str):
            miller = [int(k) for k in list(params['miller'])]

    return params

def get_one_info_from_struct_dict(struct_dict, name):

    # Simply read a dictionary key
    if isinstance(name, str):
        info = struct_dict[name]
    
    # You can have multiple nested keys
    elif isinstance(name, list):
        info = struct_dict.copy()
        for n in name:
            info = info[n]  # Read one key after the other

    else:
        ReadParamsError('Error in reading struct_dict, name is wrong.')

    return info

def get_multiple_info_from_struct_dict(struct_dict, name):

    # Extract many info at the same time 
    if isinstance(name, list):
        if all([isinstance(x, list) for x in name]):
            info = []
            for n in name:
                info.append(get_one_info_from_struct_dict(struct_dict, n))

        else:
            info = get_one_info_from_struct_dict(struct_dict, name)
    
    return info

def write_one_dict_for_db(data, name):
    """
    Prepare the dictionary to be stored in the database.
    
    """

    # Simply write a dictionary with a single key
    if isinstance(name, str):
        d = {name: data}

    # You can have multiple nested keys
    elif isinstance(name, list):
        d = {name[-1]: data}   
        for key in name[-2::-1]:
            d = {key: d}
    
    else:
        WriteParamsError('Error in writing data, name is wrong.')
    
    return d

def write_multiple_dict_for_db(data, name):
    """
    Prepare a dictionary to be stored in the database. This is a wrapper of
    write_one_dict_for_db and works with a list of data and name.

    """

    # Extract many info at the same time 
    if isinstance(name, list):
        if all([isinstance(n, list) for n in name]) and isinstance(data, list):
            d = []
            for i, n in enumerate(name):
                d.append(write_one_dict_for_db(data[i], n))

        else:
            WriteParamsError('Error in writing data, data or name is wrong.')
    
    return d


def retrieve_from_db(db_file, mp_id, collection, database=None, 
                     miller=None, name=None, is_slab=False, pymatgen_obj=True):

    # Call the navigator for retrieving structure
    nav = Navigator(db_file=db_file, high_level=database)
    
    # Define the filter to be used
    filter = {'mpid': mp_id}
    if miller is not None:
        filter.update({'miller': miller})
    
    # Extract data from the database
    structure = nav.find_data(collection=collection, filter=filter)
    
    if structure is not None:
        if name is not None:
            try:
                structure = get_one_info_from_struct_dict(structure, name)
            except:
                structure = None
        
        if pymatgen_obj and structure is not None:
            func = Slab if is_slab else Structure
            structure = func.from_dict(structure)

    return structure

def retrieve_from_tag(db_file, collection, tag, name=None, database=None):

    # Call the navigator and retrieve the simulation data from tag
    nav = Navigator(db_file=db_file, high_level=database)
    vasp_calc = nav.find_data(collection, {'task_label': tag})
    
    # Retrieve the correct dictionary and obtain the structure
    info = None
    if name is not None:
        info = get_multiple_info_from_struct_dict(vasp_calc, name)
    
    return vasp_calc, info

def orient_bulk(structure, miller, thickness, primitive=False, lll_reduce=True, 
                in_unit_planes=True):
    """
    Orient a bulk unit cell along a certain miller direction
    """
    
    # Generate the oriented bulk
    slabgen = SlabGenerator(initial_structure = structure,
                            miller_index = miller,
                            primitive=primitive,
                            lll_reduce=lll_reduce,
                            in_unit_planes=in_unit_planes,
                            min_slab_size=thickness,
                            min_vacuum_size=0)
    
    bulk_miller = slabgen.oriented_unit_cell()

    return bulk_miller

def generate_slabs(structure, miller, thickness, vacuum, thick_bulk=12,
                   center_slab=True, primitive=False, lll_reduce=True,
                   in_unit_planes=True,  ext_index=None, bonds=None, ftol=0.1, 
                   tol=0.1, repair=False, max_broken_bonds=0, symmetrize=False):
    """
    Create and return a list of slabs out of a structure.

    """

    slabs = []
    for hkl, thk, vac in zip(miller, thickness, vacuum):

        # Oriented bulk case
        if thk == 0:
            s = orient_bulk(structure, miller, thick_bulk, in_unit_planes)
        
        # Slab case
        else:
            slabgen = SlabGenerator(initial_structure=structure,
                                    miller_index=hkl,
                                    center_slab=center_slab,
                                    primitive=primitive,
                                    lll_reduce=lll_reduce,
                                    in_unit_planes=in_unit_planes,
                                    min_slab_size=thk,
                                    min_vacuum_size=vacuum)
        
            # Select the ext_index-th slab from the list of possible slabs
            s = slabgen.get_slabs(bonds=bonds, 
                                  ftol=ftol, 
                                  tol=tol, 
                                  repair=repair,
                                  max_broken_bonds=max_broken_bonds, 
                                  symmetrize=symmetrize)
            
            if ext_index:
                s = s[ext_index]

        slabs.append(s)
    
    return slabs

def _select_func_from_dict(struct_kind):

    if struct_kind == 'bulk':
        func = Structure
    elif struct_kind == 'slab':
        func = Slab
    else:
        ValueError("Wrong argument: struct_kind. Allowed values: "
                   "'bulk', 'slab'. Given value: {}".format(struct_kind)) 
    
    return func
