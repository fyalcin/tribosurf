#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:22:14 2021

Collection of very general Firetasks that can be used as constituting elements
for any advanced Workflow.

The modules contains the following Firetasks:

** DFT Simulations **:

    - FT_RelaxStructure
    General Firetask to relax a given structure, either bulk, slab or interface.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

** Database Interactions **

    - FT_MoveTagResults
    Move a subdictionary containing some results of interest from a location in
    the database to another one.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 22nd, 2021'

import os

from pprint import pprint, pformat
from pymatgen.io.vasp.inputs import Poscar
from atomate.vasp.fireworks.core import OptimizeFW, ScanOptimizeFW
from atomate.vasp.powerups import add_modify_incar
from fireworks import (
    Firework, 
    Workflow, 
    FWAction, 
    FiretaskBase, 
    FileWriteTask, 
    explicit_serialize
)
from triboflow.utils.database import Navigator
from triboflow.utils.utils import (
    read_runtask_params,
    read_default_params,
    write_multiple_dict,
    select_struct_func,
    retrieve_from_db,
    retrieve_from_tag

)
from triboflow.utils.vasp_tools import GetCustomVaspRelaxSettings
from triboflow.utils.file_manipulation import copy_output_files
from triboflow.utils.errors import RelaxStructureError, MoveTagResultsError

currentdir = os.path.dirname(__file__)


# ============================================================================
# Firetasks for DFT simulations
# ============================================================================

@explicit_serialize
class FT_RelaxStructure(FiretaskBase):
    """
    Retrieve a structure based on mp_id and entry out of a collection found in
    db_file and database. The slab is relaxed following the 'relax_type' procedure
    using an Atomate workflow. The result is stored in the same database with 
    a tag that can be provided by the user.
    
    """

    required_params = ['mp_id', 'functional', 'collection', 'entry', 'tag']
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
        structure, is_done = self.getcheck_struct(p, pymatgen_obj=False)

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
                                     entry=p['entry'],
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


# ============================================================================
# Firetasks to save, move, retrieve data from the Database
# ============================================================================

@explicit_serialize
class FT_MoveTagResults(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['mp_id', 'collection_from', 'collection_to', 'tag']
    optional_params = ['db_file', 'database_from', 'database_to', 'miller',
                       'entry_check', 'entry_to', 'entry_from', 'struct_kind', 
                       'override', 'tag_key', 'cluster_params']

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
                                default_key="MoveTagResults")

        # Check if a structure is already present in entry and check_key
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
        MoveTagResultsError.check_collection(p['collection_to'])

        if p['entry_check'] is not None:
            # Retrieve the structure from the Database
            structure = retrieve_from_db(db_file=p['db_file'], 
                                         database=p['database'], 
                                         collection=p['collection'], 
                                         mp_id=p['mp_id'],
                                         miller=p['miller'],
                                         entry=p['entry_check'],
                                         pymatgen_obj=False)
            # Check if the calculation is already done
            is_done = False if (structure is None or not bool(structure)) else False
            
        else:
            is_done = False

        return is_done
    
    def get_results_from_tag(self, p):
        
        # Retrieve the vasp_calc_output and the info
        vasp_calc, info = retrieve_from_tag(db_file=p['db_file'],
                                            collection=p['collection_from'],
                                            tag=p['tag'],
                                            tag_key=p['tag_key'],
                                            entry=p['entry_from'],
                                            database=p['database_from'])
        return vasp_calc, info
    
    def store_results(self, info, p):

        # Dictionaries are stored in entry_to[i]/entry_from[i]
        entry = []
        for n, n_t in zip(list(p['entry_to']), list(p['entry_from'])):
            entry.append(list(n).append(list(n_t)[-1]))
        
        # Prepare the list of dictionaries to be stored in the database
        info_dict = write_multiple_dict(info, entry)
    
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

        func = select_struct_func(p['struct_kind'])
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
