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
    
    - convert_dict_to_mongodb
    
    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna
    Credits: Code readapted from StackOverflow, Creative Commons licence 3.0
    https://creativecommons.org/licenses/by-sa/3.0/
    https://stackoverflow.com/questions/29267519/mongodb-update-dictionary-in-document

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

def convert_dict_to_mongodb(input_dict):
    """
    Convert a dictionary to be suitable to the update_one method of pymongo.
    This is necessary in order to avoid updating a dictionary by substituting
    the existing one with the new one.

    Parameters
    ----------
    input_dict : dict
        Dictionary to be converted to the MongoDB style to update 

    Returns
    -------
    output_dict : dict
        Output dictionary which is now compliant to MongoDB. In this way you
        can update entries containing dictionaries or nested dictionaries
        without overwriting them.
    
    Examples
    --------
    Normally when you update nested dictionaries you substitute the first one
    with the second, without having a "proper" updating, as you might intend.
    
    >>> dict1 = {'key': {'data_1': 5 }}
    >>> dict2 = {'key': {'data_2': 2 }}
    >>> dict2.update(dict1)
    >>> dict2
    {'key': {'data_1': 5}}
    
    If your intention was to have {'key': {'data_1': 5, 'data_2': 2}} the
    operation failed. The same thing occur when updating a MongoDB entry.
    To really update an existing (nested) dictionary and not overriding it run:

    >>> convert_dict_to_mongodb(dict1)
    {'key.data_1': 5}

    """
    
    output_dict = {}
    
    for key, val in input_dict.items():
        
        if not isinstance(val, dict):
            output_dict[key] = val
            continue

        for sub_key, sub_val in val.items():
            new_key = '{}.{}'.format(key, sub_key)
            output_dict[new_key] = sub_val
            if not isinstance(sub_val, dict):
                continue

            output_dict.update(convert_dict_to_mongodb(output_dict))
            if new_key in output_dict:
                del output_dict[new_key]

    return output_dict

def write_one_dict(data, entry, to_mongodb=True):
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

    to_mongodb : bool, optional
        Decide whether to return a dictionary which is conformal to the queries
        of MongoDB ($set), in order to really update the fields and not 
        overwrite dictionaries.

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
    
    # Convert the dictionary to suit MongoDB query
    if to_mongodb:
        d = convert_dict_to_mongodb(d)
    
    return d

def write_multiple_dict(data, entry, to_mongodb=True):
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
    
    to_mongodb : bool, optional
        Decide whether to return a dictionary which is conformal to the queries
        of MongoDB ($set), in order to really update the fields and not 
        overwrite dictionaries.

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
                # Convert the dictionary to suit MongoDB query
                if to_mongodb:
                    d.append(write_one_dict(data[i], n, to_mongodb))

        else:
            d = write_one_dict(data, entry, to_mongodb)

    else:
        d = write_one_dict(data, entry, to_mongodb)
    
    return d


# ============================================================================
# Retrieve structure and VASP output from DB
# ============================================================================

def retrieve_from_db(mp_id, collection, db_file=None, database=None, 
                     miller=None, entry=None, is_slab=False, pymatgen_obj=True):
    """
    [summary]

    Parameters
    ----------
    mp_id : str
        MP ID from the Materials Project, identify a material.

    collection : str
        Collection in the database to parse through.

    db_file : str or None, optional
        Location of the database. If it is None, it will be searched for a
        'localhost' on the hosting machine. The default is None.

    database : str or None, optional
        Database toquery. The default is None.

    miller : list, optional
        Miller index identifying the orientation of the slab. If it is not 
        None, it can be used as an option as filter together with `mp_id`.
        The default is None.

    entry : str or list or None, optional
        Key or list of keys to be used to extract a piece of information from 
        the `vasp_calc` dictionary. The default is None.

    is_slab : bool, optional
        Recognize the type of structure to convert a Structure or Slab
        dictionary back to a pymatgen object, by applying the `.from_dict`
        method. Meaninful only when `pymatgen_obj=True`. The default is False.

    pymatgen_obj : bool, optional
        Decide to return a pymatgen object or a dictionary. The default is 
        True.

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

def retrieve_from_tag(collection, tag, tag_key='task_label', entry=None, 
                      db_file=None, database=None):
    """
    Retrieve a dictionary field out of the database based on the combination
    {tag_key : tag} as filter. Useful to retrieve quickly the results of a 
    VASP simulation run with Atomate and typically stored in the low level DB.

    Parameters
    ----------
    collection : str
        Collection in the database to parse through.

    tag : any python object
        Object to be found within the database to identify the correct field.
    
    tag_key : str, optional
        Dict key to filter the fields of the dictionary retrieved from the 
        database. The default is 'task_label'.
        
    entry : str or list or list of lists or None, optional
        Key or list of keys to be used to extract a piece of information or 
        multiple values from the `vasp_calc` dictionary. The default is None.

    db_file : str or None
        Location of the database. If it is None, it will be searched for a
        'localhost' on the hosting machine. The default is None.

    database : str or None, optional
        Database toquery. The default is None.

    Returns
    -------
    vasp_calc : dict
        Dictionary field retrieved from `database` in `db_file`.
    
    info : any python object
        Content of the `vasp_calc` dictionary, read using entry as a key or
        a list of nested keys. If `entry=None` it returns None.

    """

    # Call the navigator and retrieve the simulation data from tag
    nav = Navigator(db_file=db_file, high_level=database)    
    vasp_calc = nav.find_data(collection, {tag_key: tag})
    
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
