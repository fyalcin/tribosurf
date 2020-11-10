import io
import os
import pymongo
from PIL import Image
from pymatgen.ext.matproj import MPRester
from atomate.vasp.database import VaspCalcDb

def GetDBJSON():
    """Return the full path to the db.json file to connect to the database.
    
    Uses the environment variable that should be set during the configuration
    of the virtual python environment for triboflow to find the config
    directory which should contain the db.json file.

    Raises
    ------
    SystemError
        If the environment variable is not set, an error will be raised.

    Returns
    -------
    db_file : str
        Full path to the db.json file that us used to access the database.

    """
    if 'FW_CONFIG_FILE' in os.environ:
        conf_file = os.environ['FW_CONFIG_FILE']
        conf_path = conf_file.rstrip('FW_config.yaml')
        db_file = conf_path + 'db.json'
        return db_file
    else:
        raise SystemError('Could not find "FW_CONFIG_FILE" environment '
                          'variable.\nPlease make sure that your python'
                          'environment is configured correctly.')

def GetDB(db_file):
    """Connect to the MongoDB database specified in the db_file.
    
    Parameters
    ----------
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database.

    Returns
    -------
    db : pymongo.database.Database
        MongoDB database accessed by pymongo

    """
    VaspCalcDB = VaspCalcDb.from_db_file(db_file)
    return VaspCalcDB.db

def ConvertBytesToImage(bytes_object):
    """Convert an image saved as bytes for storing in MongoDB to an PIL image.

    Parameters
    ----------
    bytes_object : bytes
        Image converted to bytes for storage in a MongoDB database.

    Returns
    -------
    pil_img : PIL.PngImagePlugin.PngImageFile
        Image in python image library format ready to be viewed or saved.

    """
    pil_img = Image.open(io.BytesIO(bytes_object))
    return pil_img

def ConvertImageToBytes(path_to_fig):
    """Convert an image to bytes for starage in MongoDB database.

    Parameters
    ----------
    path_to_fig : Str
        Path to the input figure.

    Returns
    -------
    bytes
        image encoded in bytes ready for storage in MongoDB.

    """
    im = Image.open(path_to_fig)
    image_bytes = io.BytesIO()
    im.save(image_bytes, format='png')
    return image_bytes.getvalue()

def GetHighLevelDB(db_file):
    """Return the triboflow MongoDB database.
    
    Parameters
    ----------
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database.
        The default is '/home/mwo/FireWorks/config/db.json'.
    
    Returns
    -------
    db : pymongo.database.Database
        Triboflow high-level MongoDB database accessed by pymongo.
    """
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    
    return tribo_db

def AddBulkToDB(structure, mp_id, db_file, functional):
    """Insert a bulk structure to the high level database.

    This will save the input structure in the database alongside some metadata
    if it is not already found there. Otherwise it will print a warning.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Structure that should be serialized and saved in the DB.
    mp_id : str
        Materials project ID number of the structure.
    db_file : str
        path to the db.json file that is used to connect to the database.
    functional : str
        Functional (PBE or SCAN) to select the correct database collection.

    Returns
    -------
    None.

    """
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.bulk_data'
    col = tribo_db[col_name]
    formula = structure.composition.reduced_formula
    
    if col.find_one({'mpid': mp_id}):
        print('{} bulk can not be added to bulk_data collection because a '
            'material with MP-ID {} is already present in the {} '
            'collection!'.format(formula, mp_id, functional))
        return
    
    col.insert_one({'mpid': mp_id,
                    'formula': formula,
                    'structure_fromMP': structure.as_dict()})
    return

def GetBulkFromDB(mp_id, db_file, functional):
    """Retrieve information about a bulk from the high-level database.

    Parameters
    ----------
    mp_id : str
        Materials project ID number of the structure to fetch.
    db_file : str
        path to the db.json file that is used to connect to the database.
    functional : str
        Functional (PBE or SCAN) to select the correct database collection.

    Raises
    ------
    IOError
        If there is no entry in the high-level db for the mp_id given, an error
        is thrown.

    Returns
    -------
    out_dict : TYPE
       Dictionary that contains all the information stored for the bulk in the
       high level database.

    """
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.bulk_data'
    col = tribo_db[col_name]
    out_dict = col.find_one({'mpid': mp_id})
    if out_dict:
        return out_dict
    else:
        raise IOError('No bulk material with MP-ID {} is found in the'
                    '{}.bulk_data collection.'.format(mp_id, functional))

def GetSlabFromDB(mp_id, db_file, miller, functional):
    """Retrieve information about a slab from the high-level database.

    Parameters
    ----------
    mp_id : str
        Materials project ID number of the structure to fetch.
    db_file : str
        path to the db.json file that is used to connect to the database.
    miller : list of int
        Miller indices for the slab as a list of three integers.
    functional : str
        Functional (PBE or SCAN) to select the correct database collection.

    Raises
    ------
    IOError
        If there is no entry in the high-level db for the mp_id-miller pair
        given, an error is thrown.

    Returns
    -------
    out_dict : TYPE
       Dictionary that contains all the information stored for the slab in the
       high level database.

    """
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.slab_data'
    col = tribo_db[col_name]
    out_dict = col.find_one({'mpid': mp_id, 'miller': miller})
    if out_dict:
        return out_dict
    else:
        raise IOError('No slab with MP-ID {} and miller indices {} was found'
                    ' in the {}.slab_data collection.'
                    .format(mp_id, miller, functional))

def GetInterfaceFromDB(name, db_file, functional):
    """Retrieve information about a slab from the high-level database.

    Parameters
    ----------
    name : str
        Unique name of the interface as determined from mp_ids and miller
        inidices by triboflow.utils.structure_manipulation.InterfaceName
    db_file : str
        path to the db.json file that is used to connect to the database.
    functional : str
        Functional (PBE or SCAN) to select the correct database collection.

    Raises
    ------
    IOError
        If there is no entry in the high-level db for the name
        given, an error is thrown.

    Returns
    -------
    out_dict : TYPE
       Dictionary that contains all the information stored for the interface
       in the high level database.

    """
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.interface_data'
    col = tribo_db[col_name]
    out_dict = col.find_one({'name': name})
    if out_dict:
        return out_dict
    else:
        raise IOError('No interface with name {} was found in the '
                    '{}.interface_data collection.'.format(name, functional))
    
def GetLastBMDatafromDB(formula, db_file='/home/mwo/FireWorks/config/db.json'):
    """Query the FireWorks database for the last bulk modulus data for formula.
    
    A query is made to the eos collection in the Fireworks database for a
    given chemical formula. The results are sorted ascending by creation date
    and than the last one is returned as a dictionary.

    Parameters
    ----------
    formula : str
        Chemical formula for which the query is made.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database.
        The default is '/home/mwo/FireWorks/config/db.json'.

    Returns
    -------
    dict
        Data generated by the get_wf_bulk_modulus workflow of atomate.

    """
    db = GetDB(db_file)
    cursor = db.eos.find({'formula_pretty': formula}).sort(
        'created_at', pymongo.DESCENDING)
    return cursor[0]

def GetPropertyFromMP(MP_ID, prop):
    """Get a certain property for a single material from the MP database.
    
    Return the selected property of the selected material from the Material
    Projects database using the MPRester from Pymatgen.

    Parameters
    ----------
    MP_ID : str
        Valid materials project ID.
    prop : str
        Property for which to query. If not in the supported list, a warning
        will be printed.

    Returns
    -------
    variable type
        The query result from the MP database. Type depends on the query.
        If nothing is found or the property is not in the list, NoneType is
        returned.
    """
    with MPRester() as mpr:
        
        if prop not in mpr.supported_task_properties:
            print('{} is not in the list of supported task properties!\n'
                'This is the supported list:\n{}'
                .format(prop, mpr.supported_task_properties))
            return None
        else:
            query=mpr.query(criteria={'material_id': MP_ID},
                            properties=[prop])
            return query[0][prop]

def GetLowEnergyStructure(chem_formula, MP_ID=None, PrintInfo=False):
    """Search MaterialsProjects for structure.
    
    A function that searches the MaterialsProject Database
    for structures that match the given chemical formula
    and selcts the one with the lowest formation energy
    per atom. If MP_ID is givem, the structure with that mp-id will
    be returned.
    
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