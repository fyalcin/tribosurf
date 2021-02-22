#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:49:54 2021

This module contains core functionalities and features common to almost all
the Firetasks. It provides tools to manage the input parameters of a Firetask
and to manipulate the dictionaries to be retrieved and stored to the database

The module contains the following functions:

** Read input parameters **:
    - read_json
    - read_runtask_params
    - read_default_params

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

** Manipulate dictionaries **:
    - get_one_info_from_struct_dict
    - get_multiple_info_from_struct_dict
    - write_one_dict_for_db
    - write_multiple_dict_for_db

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

** Interaction with Database **:
    - retrieve_from_db
    - retrieve_from_tag

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

** Input/function handling **:
    - create_tags
    - get_miller_str
    - select_struct_func

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 22nd, 2021'

import json
from uuid import uuid4

from pymatgen.core.surface import Structure, Slab
from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator
from triboflow.utils.errors import ReadParamsError, WriteParamsError


# ============================================================================
# Read input parameters from dictionaries and in Firetasks
# ============================================================================

def read_json(jsonfile):
    """
    Shortcut to easily read a json file.
        
    """
    
    with open(jsonfile, 'r') as f:
        data = json.load(f)  
    return data

def read_runtask_params(obj, fw_spec, required_params, optional_params,
                        default_file, default_key):
    """
    Read the input arguments passed by the user when initializing a Firetask
    instance. It returns a dictionary containing all the required and optional
    parameters to be used within the `run_task` method for further processing.
    Missing data is substituted with default values.

    Parameters
    ----------
    obj : Firetask object
        Instance of a Firetask-derived class, i.e. inherited from FiretaskBase.
        It is necessary to read the arguments passed as input by the user.

    fw_spec : dict
        JSON spec including all the information needed to bootstrap the job.

    required_params : dict
        List of required parameters to the Firetask.

    optional_params : dict
        List of optional parameters to the Firetask.

    default_file : str
        Path to the JSON file containing default values for `optional_params`.

    default_key : str
        Key of the JSON dict containing the default parameters of the Firetask.

    Returns
    -------
    dict
        Dictionary containing all the parameters, the optional values which have
        not been provided by the user are substituted with defaults.

    """

    defaults = read_json(default_file)
    params = {}

    # Read required and optional parameters
    for key in required_params:
        params[key] = obj[key]
    for key in optional_params:
        params[key] = obj.get(key, defaults[default_key][key])

    # Clean db_file entry
    if 'db_file' in params.keys():
        if params['db_file'] is None:
            params['db_file'] = env_chk('>>db_file<<', fw_spec)

    # Clean miller index entry
    if 'miller' in params.keys():
        if isinstance(params['miller'], str):
            params['miller'] = [int(k) for k in list(params['miller'])]

    return params

def read_default_params(default_file, default_key, dict_params):
    """
    Read the default argument read from a JSON file and compare them with the
    keys of a dictionary. If some data is missing, is substituted with default
    values. It is a generalization of `read_runtask_params`, because it takes a
    dictionary as input and can be used outside Firetasks.

    Parameters
    ----------
    default_file : str
        Path to the JSON file containing the default values.

    default_key : str
        Key of the JSON dict containing the default parameters.

    dict_params : dict
        Dictionary containing the parameters to be updated with defaults. It
        should be a subset of the dictionary extracted from `default_file`.

    Returns
    -------
    dict
        Final dictionary containing the elements of `dict_params` if present,
        else the default values read from JSON file.

    Raises
    ------
    ReadParamsError
        When unknown keys are present ind `dict_params`.
    """

    # Read the JSON file with defaults and extract the corresponding key
    defaults = read_json(default_file)
    defaults = defaults[default_key]
    params = {}

    # Check if there are some unknown parameters
    if not set(dict_params.keys()).issubset(set(defaults.keys())):
        raise ReadParamsError("The values passed as dictionary params are not "
                              "known. Allowed values for {}, read in {}: {}"
                              .format(default_key, default_file, defaults.keys()))

    # Set the parameters, missing parameters are substituted with defaults
    for key, value in defaults.items():
        params[key] = dict_params.get(key, value)
    
    return params


# ============================================================================
# Manipulate dictionary to retrieve data or create ones to be stored in DB
# ============================================================================

def get_one_info_from_dict(input_dict, entry):
    """
    Extract subdictionaries or values from an input dictionary by providing a
    list containing the key to be read one after the other.

    Parameters
    ----------
    input_dict : dict
        Input dictionary to be read.

    entry : str or list
        Contains the key(s) that should be read in sequence from the dictionary.
        The keys should be "innested" within the dictionary or an error raises.
        
        Example:
            1.  entry = 'key'
                function returns input_dict[key]

            2.  entry = ['key1', 'key2', 'key3']
                function returns input_dict['key1']['key2']['key3']

    Returns
    -------
    dict
        Subdictionary containing the desired information.
    """

    # Simply read a dictionary key
    if isinstance(entry, str):
        info = input_dict[entry]
    
    # You can have multiple nested keys
    elif isinstance(entry, list):
        info = input_dict.copy()
        for n in entry:
            info = info[n]  # Read one key after the other

    else:
        ReadParamsError('Error in reading input_dict, entry is wrong.')

    return info

def get_multiple_info_from_dict(input_dict, entry):
    """
    Extract multiple values from an input dictionary by providing a list of
    lists where each element contains the key(s) to be read.

    Parameters
    ----------
    input_dict : dict
        Input dictionary to be read.

    entry : str or list or list of lists
        Contains the key(s) that should be read in sequence from the dictionary.
        The keys should be "innested" within the dictionary or an error raises.
        If a list of lists is passed the dictionary is read multiple times.
        
        Example:
            1.  entry = 'key' or ['key1', 'key2']
                function returns:
                    input_dict[key] or input_dict['key1']['key2']

            2.  entry = [['key1', 'key2'], ['key3']]
                function returns:
                    [ 
                        input_dict['key1']['key2'],
                        input_dict['key3']
                    ]

    Returns
    -------
    dict or list of dicts
        Subdictionary containing the desired information.
    """

    # Extract many info at the same time 
    if isinstance(entry, list):
        if all([isinstance(x, list) for x in entry]):
            info = []
            for n in entry:
                info.append(get_one_info_from_dict(input_dict, n))

        else:
            info = get_one_info_from_dict(input_dict, entry)
    
    else:
        info = get_one_info_from_dict(input_dict, entry)
    
    return info

def write_one_dict(data, entry):
    """
    Prepare an object by creating a dictionary with nested keys as provided by 
    entry. Useful to create a dictionary to update a Database fields.

    Parameters
    ----------
    data : any python object
        Python object containing the data to be the value of the dictionary.

    entry : str or list
        Contains the key(s) to be used to create the dictionary. If a list is
        passed, then the elements of the list will constitute the keys for each
        inner subdictionary that is created. At innermost level is placed data.

        Example:
            1.  entry = 'key'
                function returns:
                    {'key': data}

            2.  entry = ['key1', 'key2', 'key3']
                function returns:
                    { 'key1':
                        { 'key2':
                            {
                                'key3': data
                            }
                        }
                    }

    Returns
    -------
    dict
        Dictionary containing the data.
    """

    # Simply write a dictionary with a single key
    if isinstance(entry, str):
        d = {entry: data}

    # You can have multiple nested keys
    elif isinstance(entry, list):
        d = {entry[-1]: data}   
        for key in entry[-2::-1]:
            d = {key: d}
    
    else:
        WriteParamsError('Error in writing data, entry is wrong.')
    
    return d

def write_multiple_dict(data, entry):
    """
    Prepare multiple dictionaries containing the passed data. It is a wrapper of
    `write_one_dict` and works with a list of data as entry.

    Parameters
    ----------
    data : list of any python object
        List of Python object containing to be the value of the dictionaries.

    entry : str or list of lists
        Contains the keys to be used to create the dictionaries. If a list of
        lists is passed, then multiple dictionaries are created and returned as
        a list. In that case: `len(data) == len(entry)`, otherwise data will be
        considered as a list object and placed as a single dictionary element.

        Example:
            1.  entry = 'key'
                function returns:
                    {'key': data}

            2.  entry = [['key1', 'key2'], ['key3']]
                function returns:
                [
                    { 'key1':
                        { 
                            'key2': data[0]
                        }
                    },

                    {
                        'key3': data[1]
                    }
                ]

    Returns
    -------
    dict or list of dicts
        Dictionaries containing the data.
    """

    # Extract many info at the same time 
    if isinstance(entry, list):
        if all([isinstance(n, list) for n in entry]) and isinstance(data, list) and len(data) == len(entry):
            d = []
            for i, n in enumerate(entry):
                d.append(write_one_dict(data[i], n))

        else:
            d = write_one_dict(data)

    else:
        d = write_one_dict(data)
    
    return d


# ============================================================================
# Retrieve structure and VASP output from DB
# ============================================================================

def retrieve_from_db(db_file, mp_id, collection, database=None, 
                     miller=None, entry=None, is_slab=False, pymatgen_obj=True):
    """
    [summary]

    Parameters
    ----------
    db_file : [type]
        [description]
    mp_id : [type]
        [description]
    collection : [type]
        [description]
    database : [type], optional
        [description], by default None
    miller : [type], optional
        [description], by default None
    entry : [type], optional
        [description], by default None
    is_slab : bool, optional
        [description], by default False
    pymatgen_obj : bool, optional
        [description], by default True

    Returns
    -------
    [type]
        [description]
    """

    # Call the navigator for retrieving structure
    nav = Navigator(db_file=db_file, high_level=database)
    
    # Define the filter to be used
    filter = {'mpid': mp_id}
    if miller is not None:
        filter.update({'miller': miller})
    
    # Extract data from the database
    structure = nav.find_data(collection=collection, filter=filter)
    
    if structure is not None:
        if entry is not None:
            try:
                structure = get_one_info_from_dict(structure, entry)
            except:
                structure = None
        
        if pymatgen_obj and structure is not None:
            func = Slab if is_slab else Structure
            structure = func.from_dict(structure)

    return structure

def retrieve_from_tag(db_file, collection, tag, entry=None, database=None):
    """
    [summary]

    Parameters
    ----------
    db_file : [type]
        [description]
    collection : [type]
        [description]
    tag : [type]
        [description]
    entry : [type], optional
        [description], by default None
    database : [type], optional
        [description], by default None

    Returns
    -------
    [type]
        [description]
    """

    # Call the navigator and retrieve the simulation data from tag
    nav = Navigator(db_file=db_file, high_level=database)
    vasp_calc = nav.find_data(collection, {'task_label': tag})
    
    # Retrieve the correct dictionary and obtain the structure
    info = None
    if entry is not None:
        info = get_multiple_info_from_dict(vasp_calc, entry)
    
    return vasp_calc, info


# ============================================================================
# Secondary tools to work with input parameters and functions
# ============================================================================

def create_tags(prefix):
    """
    Create a tag out of a prefix.

    """

    # Create a list of tags
    if isinstance(prefix, list):
        tag = [n + '_' + str(uuid4()) for n in prefix]

    else:
        tag = prefix + '_' + str(uuid4())
    
    return tag

def get_miller_str(miller):
    """
    Convert a miller index from list to string.

    """
    return ''.join(str(s) for s in miller)

def select_struct_func(struct_kind):
    """
    Select a function to work on pymatgen structure, depending on the value
    of `struct_kind` it returns either `pymatgen.core.surface.Structure` or 
    `pymatgen.core.surface.Slab`.

    """

    if struct_kind == 'bulk':
        func = Structure
    elif struct_kind == 'slab':
        func = Slab
    else:
        ValueError("Wrong argument: struct_kind. Allowed values: "
                   "'bulk', 'slab'. Given value: {}".format(struct_kind)) 
    return func
