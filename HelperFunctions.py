"""A collection of HelperFunctions to be used for the FireFlow project."""


def GetCustomVaspRelaxSettings(comp_parameters, relax_type):
    """Make custom vasp settings for relaxations.
    

    Parameters
    ----------
    comp_parameters : dict
        Computational parameters dictionary which is usually created partly
        by the user and then filled automatically with defaults and after
        convergence tests.
    relax_type : str
        Specifies what is to be relaxed in what way. Check 'allowed_types'
        for a list of choices.

    Raises
    ------
    SystemExit
        If a non-supported relax_type is passed, the process terminates.

    Returns
    -------
    vis : str
        A vasp input set for pymatgen.
    uis : dict
        User input settings that will override the standard setting in the vis.
    vdw : str
        Specifies which vdw functional is used. (optB86b or rVV10)

    """
    allowed_types = ['bulk_full_relax', 'bulk_vol_relax', 'bulk_pos_relax',
                     'bulk_shape_relax',
                     'slab_shape_relax', 'slab_pos_relax',
                     'interface_shape_relax', 'interface_pos_relax']
    
    if relax_type not in allowed_types:
        raise SystemExit('relax type is not known. Please select from: {}'
                         .format(allowed_types))
    
    #Set user incar settings:
    uis = {}
    uis['NEDOS'] = 3001
    uis['PREC'] = 'Accurate'
    uis['GGA_COMPAT'] = '.FALSE.'
    uis['LASPH'] = '.TRUE.'
    uis['LORBIT'] = 11
    uis['MAXMIX'] = 100
    uis['NELMIN'] = 4
    
    if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
        uis['NELMDL'] = -15
        uis['EDIFFG'] = -0.05
    else:
        uis['EDIFFG'] = -0.01
    
    if relax_type.endswith('full_relax'):
        uis['ISIF'] = 3
    elif relax_type.endswith('pos_relax'):
        uis['ISIF'] = 2
    elif relax_type.endswith('vol_relax'):
        uis['ISIF'] = 7
    elif relax_type.endswith('shape_relax'):
        uis['ISIF'] = 5
    
    if 'encut' in comp_parameters:
        uis['ENCUT'] = comp_parameters['encut']
        
    if 'use_spin' in comp_parameters:
        if comp_parameters['use_spin']:
            uis['ISPIN'] = 2
        else:
            uis['ISPIN'] = 1
    
    if 'is_metal' in comp_parameters:
        if comp_parameters['is_metal']:
            uis['SIGMA'] = 0.2
            uis['ISMEAR'] = 2
        else:
            uis['SIGMA'] = 0.05
            uis['ISMEAR'] = -5
    else:
        uis['SIGMA'] = 0.1
        uis['ISMEAR'] = 0
        
        
    #set van der Waals functional. Note that as of now, 'functional' must be
    #specified for vdw to work!
    if set(('use_vdw', 'functional')) <= comp_parameters.keys():
        if comp_parameters['use_vdw']:
            if comp_parameters['functional'] == 'SCAN':
                vdw = 'rVV10'
            else:
                vdw = 'optB86b'
        else:
            vdw = None
    else:
        vdw = None
        
    if 'functional' in comp_parameters:
        if comp_parameters['functional'] == 'SCAN':
            vis = 'MPScanRelaxSet'
            #Algo All does not play well with tetrahedron method
            if 'is_metal' in comp_parameters:
                if not comp_parameters['is_metal']:
                    uis['SIGMA'] = 0.1
                    uis['ISMEAR'] = 0
        else:
            vis = 'MPRelaxSet'
    else:
        vis = 'MPRelaxSet'
        
    return vis, uis, vdw

def UpdateNestedDict(d, u):
    """Update the dictionary d with the update u without losing data in d.
    
    This function allows to update a nested dictionary on an arbitrary level
    without overwriting any data present in the key that is being updated.
    This entire function is copied from the follwing stackoverflow post:
    https://stackoverflow.com/a/3233356

    Parameters
    ----------
    d : dict
        (Nested) dictionary to be updated. 
    u : dict
        (Nested) update to the dictionary.

    Returns
    -------
    d : dict
        Updated dictionary

    """
    import collections.abc
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = UpdateNestedDict(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def WriteNestedDictFromList(key_list, value, d={}):
    """Write a value to a nested dict defined by a list of strings.
    
    A list of keys is given and a nested dict is returned. The input value
    will be put into the last key of the key_list. E.g.:
        value = 'x' and key_list = ['a', 'b', 'c'] will return
        {'a': {'b': {'c': 'x'}}}
    
    Parameters
    ----------
    key_list : list of str
        A list of keys with wich to form the nested dictionary.
    value :
        The object to put at the end of the nested dictionary Can be any type.
    d : dict, optional
        A dictionary can be passed, but it will be overwritten. It is needed
        for the recursive definition. The default is {}.

    Returns
    -------
    d : dict
        A nested dictionary with the input value at the end of it.

    """
    if len(key_list) > 1:
        d[key_list[0]] = {}
        d[key_list[0]] = d.get(key_list[0], {})
        WriteNestedDictFromList(key_list[1:], value, d[key_list[0]])
    else:
        d[key_list[0]] = value
    return d

def GetValueFromNestedDict(dictionary, key_list):
    """Take a list of keys and a nested dictionary and return the value.
    
    A (usually nested) dictionary is given along side a list of keys to this
    dictionary. The value at the end of this key_list is returned using a
    recursive method. E.g. key_list=['key_a', 'key_b', 'key_c'] will return
    dictionary['key_a']['key_b']['key_c']
    Returns <None> when there is a KeyError
    

    Parameters
    ----------
    dictionary : dict
        Input dictionary which can be nested
    key_list : list of str
        The list of keys to the nested dictionary at the end of which the
        value that is needed is located

    Returns
    -------
        Returns whatever is stored at the appropriate position in the
        dictionary. Type can vary!
        If there is a KeyError on some stage, None is returned

    """
    if len(key_list) > 1:
        if key_list[0] in dictionary:
            return GetValueFromNestedDict(dictionary[key_list[0]],
                                          key_list[1:])
        else:
            return None
    return dictionary.get(key_list[0])
        

def RemoveMatchingFiles(list_of_patterns):
    """
    Remove all files matching the patterns (wildcards) in the current
    directory.
    """
    import os, glob
    remove_list=[]
    for pattern in list_of_patterns:
        remove_list.extend(glob.glob(pattern))
    if remove_list is []:
        return
    else:
        for file_to_remove in remove_list:
            os.remove(file_to_remove)
        return

def WriteFileFromDict(Dict, Filename):
    """
    Takes an input dictionary 'Dict' and an output filename (or
    path) 'Filename' as a string and writes the dictionary as
    the output file. If a list is found as a value of the dictionary,
    its entries are printed without colons or brakets.
    E.g. {'ENCUT': 320, 'ALGO: 'FAST', 'MAGMOM': [3.0, -3.0]}
    is written as:
                ENCUT = 320
                ALGO = Fast
                MAGMOM = 3.0 -3.0
    to the file.
    """
    Out_file = []
    for key in Dict.keys():
        if type(Dict[key]) is list:
            str_list = [str(x) for x in Dict[key]]
            value = ' '.join(str_list)
        else:
            value = str(Dict[key])
        Out_file.append(str(key)+' = '+value)
    with open(Filename, 'w') as out:
        for line in Out_file:
            out.write(line+'\n')
    return

def GetGapFromMP(MP_ID):
    """Get the bandgap of a structure from the MaterialsProject database.
    
    Parameters
    ----------
    MP_ID : str
        Materials Project material_id

    Returns
    -------
    float
        Bandgap of the material (is usually not accurate, but can be used
        to decide if the material is a metal or not).

    """
    from pymatgen import MPRester
    with MPRester() as mpr:
        band_gap = mpr.query(criteria={'material_id': MP_ID},
                             properties=['band_gap'])
    return band_gap[0]['band_gap']

def GetLowEnergyStructure(chem_formula, MP_ID=None, PrintInfo=False):
    """
    A function that searches the MaterialsProject Database
    for structures that match the given chemical formula
    and selcts the one with the lowest formation energy
    per atom. If several 
    Inputs: 
    chem_formula (str): Required input of a chemical formula
                        e.g.: NaCl, Fe2O3, SiO, FeCW
    MP_ID (str):        Optional Input of Materials Project ID
                        of the exact desired structure
                        e.g. 'mp-990448'
    PrintInfo (bool):   Optional variable defaults to "False".
                        if set to "True", will print some
                        information about how many structures
                        have been found and the ID of the selcted
                        one.
                        
    Returns: pymatgen Structure object and associated MP_ID
    """
    
    from pymatgen import Structure, MPRester

    # load structure from MaterialsProjct Website
    with MPRester() as mpr:
        if MP_ID:
            struct = mpr.get_structure_by_material_id(MP_ID)
            return struct, MP_ID
        else:
            id_list = mpr.query(criteria={'pretty_formula': chem_formula,
                                        'e_above_hull': 0.0},
                              properties=['material_id'])
            if id_list == []:
                raise NameError('{} has not been found in the MaterialsProject'
                            'database'.format(chem_formula))
            else:
                MP_ID = id_list[0]['material_id']
                struct = mpr.get_structure_by_material_id(MP_ID)
                return struct, MP_ID
