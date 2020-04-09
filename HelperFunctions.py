"""A collection of HelperFunctions to be used for the FireFlow project."""


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
