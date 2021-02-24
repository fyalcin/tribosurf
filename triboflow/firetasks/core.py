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
    convert_dict_to_mongodb,
    write_multiple_dict,
    select_struct_func,
    retrieve_from_db,
    retrieve_from_tag

)
from triboflow.utils.vasp_tools import get_custom_vasp_static_settings
from triboflow.utils.file_manipulation import copy_output_files
from triboflow.utils.errors import RelaxStructureError, MoveTagResultsError

currentdir = os.path.dirname(__file__)


# ============================================================================
# Firetasks for DFT simulations
# ============================================================================

@explicit_serialize
class FT_RelaxStructure(FiretaskBase):
    """
    Retrieve a structure based on `mp_id` (and `miller`) from a location in the
    database, identified by: `db_file`, `database`, `collection`, `entry`. 
    The structure is converged setting up a VASP calculation based on Atomate
    workflows. The result is stored in the same database location and can be 
    identified by a tag object, provided by the user, associated to a 
    'task_label' entry. The type of the simulation is identified by 'relax_type'. 

    Parameters
    ----------
    mp_id : str
        MP-ID of the structure from the MP database, used to identify the
        structures in local databases.
        
    functional : str
        Functional for the pseudopotential to be adopted.

    collection : str
        Collection where the structure is present. The default collections 
        used in the "tribchem" database are identified as functional+string, 
        where string is usually: '.bulk_data', 'slab_data', 'interface_name'.

    entry : str or list
        Location of the structure to be retrieved from database and collection.

    tag : str (any python object)
        Defined by user, it is automatically associated from Atomate workflows 
        to an entry, with key 'task_label', in the first layer of the db field 
        containing the output of the VASP simulations. The results are stored in
        the Atomate database. Tag can be used to keep track of where the DFT 
        calculation results have been stored in the database. In principle it 
        can be any python object.

    db_file : str or None, optional
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.
    
    database : str, optional
        Name of the database where the structure will be retrieved. 
        The default is "tribchem".

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to `get_custom_vasp_relax_settings`. The default is 'slab_pos_relax'.

    comp_params : dict, optional
        Computational parameters to simulate the structure with DFT. If an empty 
        dictionary is passed, defaults are used. The default is {}.

    miller : list of int, optional
        Miller indexes (h, k, l) used to identify unequivocally a slab within
        the database. If nothing is passed, the material will be only identified
        by means of `mp_id`. The default is None.
    
    check_key : str or list, optional
        Check if a given data is already present within a database location, to
        understand if the DFT simulation has been already done. The location is
        identified by db_file, database, collection, entry. Once the data 
        corresponding to entry is extracted, if data['check_key'] does 
        exist then the DFT simulation is not done. If it is None, the simulation
        will be always started. The default is None.

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
        """
        Retrieve a structure from the database and eventually check if a DFT
        simulation has been already done

        Parameters
        ----------
        p : dict
            All the input parameters of the Firetasks, placed in a dictionary.

        pymatgen_obj : bool, optional
            Decide whether to directly convert the dictionary extracted to the
            database to a structure or not. The default is False.

        Returns
        -------
        structure : dict or pymatgen structure
            Dictionary containing the structure.
        
        is_done : bool
            If a simulation output is already present in the database. If it
            is True, the calculation will not be performed.

        """
        
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
        """
        Set necessary VASP parameters and create an Atomate workflow to start a 
        detour and perform a DFT simulation on the given structure structure.

        Parameters
        ----------
        structure : pymatgen structure
            Crystalline atomic structure to perform a Kohn-Sham minimization.

        p : dict
            All the input parameters of the Firetasks, placed in a dictionary.

        dfl : str
            Path to location of defaults, set missing computational parameters.

        Returns
        -------
        wf : Fireworks.Workflow
            Workflow to start a VASP simulation on your structure.
        """

        # Check tag and computational parameters
        tag = p['tag']
        comp_params = p['comp_params']
        if not bool(comp_params):
            comp_params = read_default_params(dfl, 'comp_params', comp_params)
        
        # Set options for vasp
        vis = get_custom_vasp_static_settings(structure, comp_params, p['relax_type'])
        
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
    Firetask to move some generic data from a certain location of the database 
    to another one. In principle this Firetask could be used in a very general
    way. The most common usage is to locate the results of an Atomate Workflow
    in the low level database and move output data of interest to the high level
    database of interest. 
    
    The Firetask does the following steps:

    1 - If `check_key` is not None, a control is done to understand if data is 
    already present and eventually stop the transfer to avoid overriding.
    The location of the data in the db is identified by means of: `db_file`, 
    `database_to`, `collection_to`. The MongoDB field is identified with `mp_id`
    (and `miller`) and data is retrieved. If data_dict['check_key']: stop.

    2 - The field containing the data of interest is extracted. The query is 
    done on: `db_file`, `database_from`, `collection_from`. The correct field
    is identified using the filter:
        {
            `tag_key`: `tag`
        }

    3 - Once the field dictionary has been retrieved, the data to be transferred
    is identied by `entry_from`, which could be: str, list, list of lists.
    To transfer more data it is necessary to use lists of lists, where each
    list contains the nested keys within the extracted dictionary containing 
    the data that I want to transfer.

    4 - The selected data is stored in the destination database, the location is
    identified by: `db_file`, `database_to`, `collection_to`, `entry_to`.
    `entry_to` can be again, str, list or list of lists. To identify a single 
    element it is necessary to use str or list, the latter in the case of a
    data value nested within the source or destination dictionary. However, if
    you want to transfer more data, it is mandatory to use list of lists, and
    in that case len(`entry_from`) == len(`entry_to`).
    Warning: The Firetask does not cancel the source data!

    Examples
    --------
    Here is a simple example explaining how MoveTagResults logically works.

    - Information concerning the source data and location:
        db_file = None
        database_from = "FireWorks"
        collection_from = "coll.tasks"
        entry_from = [['test', 'energy'], ['test', 'energy2']]

    - Information concerning the source tag:
        tag = True
        tag_key = "transfer_test"

    - Information concerning the destination location:
        mp_id = 'mp-126'
        db_file = None
        database_to = "tribchem"
        collection_to = "PBE.slab_data"
        entry_to = [['energy'], ['data_back', 'energy2']]
        check_key = 'is_done'

    a) The source field would be something like that:
    {
        "_id" : ...,
        "transfer_test" : true,
        "test" : {
            "energy" : 10,
            "energy2" : 50
        }
    }

    b) The destination field would be something like that:
     {
        "_id" : ...,
        "mp_id" : "mp-126"
    }

    Running the Firetasks you would:

    1. Check b-field to see if a key named 'check_key' is present. It is not, so
       the process will continue.
    2. Identify univocally the a-field to be the correct source containing the 7
       data of interest. This is done matching tag_key and key with the entry
       of the field dictionary: 'transfer_test'.
    3. Extract both energy and energy2 values from the a-field.
    4. Find the exact destination location with a filter based on mp_id and 
       place there, following the path provided by 'entry_to'.

    In the end, the b-field becomes:
     {
        "_id" : ...,
        "mp_id": "mp-126",
        "energy" : 10,
        "data_back" : {
            "energy2" : 50,
        }
    }   

    """

    required_params = ['mp_id', 'collection_from', 'collection_to', 'tag']
    optional_params = ['db_file', 'database_from', 'database_to', 'miller',
                       'check_entry', 'entry_to', 'entry_from', 'struct_kind', 
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
            wf = self.user_output(vasp_calc, p, dfl)
            if wf is not None:
                return FWAction(detours=wf, update_spec=fw_spec)
        
        else:  
            return FWAction(update_spec=fw_spec)
    
    def check_struct(self, p):
        """
        Check if there exists an entry called 'check_entry' in the destination
        field. If it is found the transfer is stopped.

        Parameters
        ----------
        p : dict
            All the input parameters of the Firetasks, placed in a dictionary.

        Returns
        -------
        is_done : bool
            If a simulation output is already present in the database. If it
            is True, the calculation will not be performed.

        """
        
        # Check if collection does exist
        MoveTagResultsError.check_collection(p['collection_to'])
        
        if p['check_entry'] is not None:
            # Retrieve the structure from the Database
            check_dict = retrieve_from_db(db_file=p['db_file'], 
                                          database=p['database_to'], 
                                          collection=p['collection_to'], 
                                          mp_id=p['mp_id'],
                                          miller=p['miller'],
                                          entry=p['check_entry'],
                                          pymatgen_obj=False)

            # Check if the calculation is already done
            is_done = True
            if check_dict is None or check_dict == {}:
                is_done = False
                
        else:
            is_done = False

        return is_done
    
    def get_results_from_tag(self, p):
        """
        Identify the correct field in `collection_from` by using tags and extract
        the output data of interests with `entry_from`.

        """
        
        # Retrieve the vasp_calc_output and the info
        vasp_calc, info = retrieve_from_tag(db_file=p['db_file'],
                                            collection=p['collection_from'],
                                            tag=p['tag'],
                                            tag_key=p['tag_key'],
                                            entry=p['entry_from'],
                                            database=p['database_from'])
        return vasp_calc, info
    
    def store_results(self, info, p):
        """
        Store the results to the destination database and collection.
        
        """
        
        # Prepare the list of dictionaries to be stored in the database
        info_dict = write_multiple_dict(info, p['entry_to'])  

        if not isinstance(info_dict, list):
            info_dict = [(info_dict)]
    
        # Prepare the database and options where to store data
        nav = Navigator(db_file=p['db_file'], high_level=p['database_to'])
        filter = {'mpid': p['mp_id']}
        if p['miller'] is not None:
            filter.update({'miller': p['miller']})        

        # Finally store the data
        for d in info_dict:
            nav.update_data(p['collection_to'], filter, {'$set': d}, upsert=True)
    
    def user_output(self, vasp_calc, p, dfl):

        # Get cluster params and set missing values
        cluster_params = p['cluster_params']
        cluster_params = read_default_params(default_file=dfl, 
                                             default_key="cluster_params", 
                                             dict_params=cluster_params)
        
        # Handle structure file output:
        # BE CAREFUL: it works only with structures elements
        if cluster_params['file_output']:
            
            # Recover the structure from the dictionary
            func = select_struct_func(p['struct_kind'])
            structure = func.from_dict(vasp_calc['output']['structure'])
            
            # Output to screen
            print('')
            print('Relaxed output structure as pymatgen.surface.Slab dictionary:')
            pprint(structure.as_dict())
            print('')
            
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
