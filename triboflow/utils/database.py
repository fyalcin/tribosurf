#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:15:11 2021

Classes to manage data from local and online DataBases at a high level.

@author: omarchehaimi
"""

__author__ = 'Omar Chehaimi'
__credits__ = 'This module is based on the Triboflow package, Michael Wolloch'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 1st, 2021'

import io
import os
from datetime import datetime
from pathlib import Path, PurePosixPath
import pymongo
from PIL import Image
from pymatgen.ext.matproj import MPRester
from atomate.vasp.database import VaspCalcDb

from triboflow.core.logging import LoggingBase
from triboflow.utils.errors import NavigatorError, StuctNavigatorError, NavigatorMPError

# Logging
configurations = LoggingBase.get_config()
logging_path = LoggingBase.logging_path()
log = LoggingBase(
    name='database',
    console_level=configurations['logging']['database_logging_level'],
    path=logging_path+'/database.log',
    file_level=configurations['logging']['database_logging_level'])

# ============================================================================
# Navigator Classes
# ============================================================================

class Navigator:
    """
    The Navigator class is a high level interface that deals with database's 
    operations. 
    In this class the following methods imported by MongoDB database are used:
        - find_one
        - find_many
        - insert_one
        - insert_many
        - update_one
        - update_many
        - delete_one
        - delete_many
        - drop_data

    PyMongo official documentation can be found here:
    https://pymongo.readthedocs.io/en/stable/api/index.html

    Attributes
    ----------
    db_file : str
        Path to the database.

    high_level : str
        High level database name.

    db : VaspCalcDb
        A VASP database type.

    path : str
        The string containing the path to the database.

    Methods
    -------
    __get_db(db_file)
        Create an instance of a VaspCalcDb database.

    __initialize_obj_collection(collection)
        Initialize a MongoDB object collection.
    
    update_data(collection, filter, new_values, )
        Update a single document matching the filter with the new value.
    
    update_many_data(collection, filter, new_values, upsert)
        Update many documents that match the filter with the new values.
    
    insert_data(collection, data, duplicates, message)
        Insert a single document in the collection.
    
    insert_many_data(collection, data, duplicates)
        Insert an iterable of documents in the collection.

    find_data(collection, filter)
        Get a single document in the collection that matches the filter.
    
    find_many_data(collection, filter)
       Get all the documents in the collection that match the filter.
    
    delete_data(collection, filter)
        Delate a single document that matches the filter.
    
    delete_many_data(collection, filter)
        Delate one or more documents that match the filter.
    
    drop_data(Drop all the documents in the collection.)
        Drop all the documents in the collection.

    """

    def __init__(self, db_file='localhost', high_level=None):
        """
        Parameters
        ----------
        db_file : str, optional
            Location where the database is saved. The default is 'localhost'.

        high_level : str or None, optional
            Decide whether to use an high level database or not, to store the
            data of the simulations. The name of that DB can be passed as a
            string to high_level. The default is None.

        """

        db, db_path = self.__get_db(db_file)
        self.db = db
        self.path = db_path
        
        if high_level is not None:
            self.db = self.db.client[high_level]

    def __get_db(self, db_file):
        """ 
        Connect to the MongoDB database specified in the db_file. 
        
        Parameters
        ----------
        db_file : str
            Path where the database file is saved.
        
        Returns
        -------
        vasp_db.db : VaspCalcDb
            VASP database object.

        db_file : str
            Database file path.

        """

        if db_file is None or db_file == 'localhost' or db_file == 'local':
            if 'FW_CONFIG_FILE' in os.environ:
                conf_file = os.environ['FW_CONFIG_FILE']
                conf_path = conf_file.rstrip('FW_config.yaml')
                db_file = conf_path + 'db.json'
            else:
                raise NavigatorError('Could not find "FW_CONFIG_FILE" environment '
                                     'variable.\nPlease make sure that your python'
                                     'environment is configured correctly.')

        try:
            vasp_db = VaspCalcDb.from_db_file(db_file)
        except: 
            raise NavigatorError('The database file does not exist in path: {}'
                                 .format(db_file))

        log.info('Successfully connected to: {}.'.format(db_file))
        return vasp_db.db, db_file
    
    def __initialize_obj_collection(self, collection):
        """
        Initialize a MongoDB object collection.

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database o r in the VaspCalcDb object.
        
        Returns
        -------
        collection_obj : VaspCalcDb
            Database object containing the collection of the database.     
        """

        if isinstance(collection, str):
            collection_obj = self.db.coll[collection]
        elif isinstance(collection, VaspCalcDb):
            collection_obj = collection
        else:
            raise NavigatorError('{} is not a valid data type. The collection '
                                 'must be a string or VaspCalcDb.'
                                 ' '.format(type(collection)))
        
        return collection_obj

    def update_data(self, collection, filter, new_values, upsert=False):
        """
        Update a single document matching the filter with the new value.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.update_one

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database or in the VaspCalcDb object.

        filter : dict
            Dictionary containing the document to update.

        new_values : dict
            The new values to update in the document.
            This dictionary must respect a particular sintax, e.g.:
                {'$set': {dictionary with the new data}}
            The '$set' variable is used for updating all the value in the 
            dictionary. Use '$inc' for integer.
            For more options check the PyMongo documentation.
        
        upsert : bool
            PyMongo parameter for the update_one function. If True update_one 
            performs an insertion if no documents match the filter.

        """

        collection_obj = self.__initialize_obj_collection(collection)

        log.info('Updating the collection {} withe the new data {}.'
                 ''.format(collection, new_values))
        collection_obj.update_one(filter, new_values, upsert)
    
    def update_many_data(self, collection, filter, new_values, upsert=False):
        """
        Update many documents that match the filter with the new values.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.update_many

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database or in the VaspCalcDb object.
        
        filter : dict
            Dictionary containing the documents to update.

        new_values : dict
            The new values to update in the document.
        
        upsert : bool
            PyMongo parameter for the update_one function. If True update_one 
            performs an insertion if no documents match the filter.

        """

        collection_obj = self.__initialize_obj_collection(collection)

        log.info('Updating the collection {} withe the new data {}.'
                 ''.format(collection, new_values))
        collection_obj.update_many(filter, new_values, upsert)

    def insert_data(self, collection, data, duplicates=False, message=None):
        """
        Insert a single document in the collection.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.insert_one

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database or a VaspCalcDb object.

        data : dict
            Data to add in the database.
        
        duplicates : bool, Optional
            If True the data are saved in the database even if they are 
            duplicates.
        
        message : str, Optional
            Custom message to write in the console and/or in the log file.

        """

        collection_obj = self.__initialize_obj_collection(collection)
        
        if not duplicates:
            if self.find_data(collection, data):
                if message:
                    log.warning(message)
                else:
                    log.warning('{} already exist in {} collection.'
                                'Use the duplicates flag for saving multiple '
                                'times the same document.'
                                ''.format(data, collection))
                return

        log.info('Writing {} in the collection {}.'.format(data, collection))
        collection_obj.insert_one(data)

    def insert_many_data(self, collection, data, duplicates=False):
        """
        Insert an iterable of documents in the collection.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.insert_many

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database or a VaspCalcDb object.

        data : list
            Many sets of data to add in the database. In the list there are all
            the dictionaries to be added in the database.

        duplicates : bool, Optional
            If True the data are saved in the database even if they are 
            duplicates.

        """

        collection_obj = self.__initialize_obj_collection(collection)

        if not duplicates:
            document_index = 0
            duplicates = []
            for document_index in range(0, len(data)):
                if self.find_many_data(collection, data[document_index]):
                    log.warning('{} already exist in {} collection.'
                                ''.format(data[document_index], collection))
                    # Save the index of all duplicate documents
                    duplicates.append(document_index)

            # Removing all the data already present in the database.
            # Starting from the last index to not change the list  
            for document_index in sorted(duplicates, reverse=True):
                data.pop(document_index)

            if len(data) == 0:
                log.warning('All the documents already exist in collection {}.'
                            ' Use the duplicates flag for saving multiple times'
                            ' the same document.'.format(collection))
                return

        log.info('Writing {} in the collection {}.'.format(data, collection))
        collection_obj.insert_many(data)

    def find_data(self, collection, filter):
        """
        Get a single document in the collection that matches the filter.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.find_one

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database or VaspCalcDb object.

        filter : dict, optional
            Dictionary containing the name of the document to find.
            If filter is empty all the documents in the collection are returned.

        Returns
        -------
        data : variable object
            A variable object depending on the type of the retrived data.
            data is None if nothing has been found.

        """

        collection_obj = self.__initialize_obj_collection(collection)

        data = collection_obj.find_one(filter)

        if not data:
            log.warning('There are no data for {}.'.format(filter))
            return data

        log.info('{} has been found in {}.'. format(filter, collection))
        return data

    def find_many_data(self, collection, filter):
        """
        Get all the documents in the collection that match the filter.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.find

        Parameters
        ----------

        collection : str, VaspCalcDb
            Name of the collection in the database or VaspCalcDb object. 

        filter : dict, optional
            Dictionary containing the name of the document to find.
            If filter is empty all the documents in the collection are returned.

        Returns
        -------
        data : variable object
            A variable object depending on the type of the data.
            data is None if nothing has been found.

        """

        collection_obj = self.__initialize_obj_collection(collection)

        data = collection_obj.find(filter)

        if not data:
            log.warning('There are no data for {}.'.format(filter))
            return data

        log.info('{} has been found in {}.'. format(filter, collection))
        return data

    def delete_data(self, collection, filter):
        """
        Delate a single document that matches the filter.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.delete_one

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database or VaspCalcDb object.

        filter : dict
            Document to be removed.

        """

        collection_obj = self.__initialize_obj_collection(collection)

        log.info('Deleting {} from the collection {}.'
                 ''.format(filter, collection))
        collection_obj.delete_one(filter)

    def delete_many_data(self, collection, filter):
        """
        Delate one or more documents that match the filter.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/collection.html#pymongo.collection.Collection.delete_many

        Parameters
        ----------
        collection : str, VaspCalcDb
            Name of the collection in the database or VaspCalcDb object.

        filter : dict
            Documents to be removed.

        """

        collection_obj = self.__initialize_obj_collection(collection)

        collection_obj.delete_many(filter)

    def drop_data(self, collection):
        """
        Drop all the documents in the collection.

        https://pymongo.readthedocs.io/en/stable/api/pymongo/database.html#pymongo.database.Database.drop_collection

        Parameters
        ----------
        collection : str
            Name of the collection to be removed.

        """
        log.critical('This will drop all entries {} from the database. '
                     'Write the current date (YYYY-mm-dd, e.g. 2021-02-01) '
                     'to confirm: '.format(filter))
        user_date = input()
        current_date = datetime.today().strftime("%Y-%m-%d")

        if user_date == current_date:
            log.critical('Removing {} from the database.'
                         ''.format(collection))
            collection_obj = self.__initialize_obj_collection(collection)
            self.db.collection_obj.drop()
        else:
            log.critical('The current date is wrong!!! '
                         'No entries in the database have been removed.')


class StructureNavigator(Navigator):
    """
    Child class of Navigator in which are implemented all the methods required 
    for writing and loading data from the database about bulk, slab, 
    and interface.
    
    These functions are taken from the TriboFlow package, Michael Wolloch.

    Attributes
    ----------
    high_level : str
        High level database name.
    
    Methods
    -------
    add_bulk_to_db(structure, mp_id, functional, message)
        Insert a bulk structure in the triboflow high level database.

    get_bulk_from_db(mp_id, functional)
        Get the data about bulk from the triboflow high level database.

    get_slab_from_db(mp_id, functional, miller)
        Get the data about slab from the triboflow high level database.
    
    get_interface_from_db(name, functional)
        Get the data about intreface from the eos database.

    get_last_bmd_data_from_db(formula)
        Get the last bulkmodule from the FireWorks database.

    """

    def __init__(self, db_file, high_level):
        super().__init__(db_file=db_file, high_level=high_level)

    def add_bulk_to_db(self, structure, mp_id, functional, message=None):
        """
        Insert a bulk structure in the triboflow high level database.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            Structure of the bulk in the pymatgen Structure format.

        mp_id : str
            Materials Project id of the structure.
        
        functional : str
            Functional. It could be PBE por SCAN.
        
        message : str, Optional
            Custom message to write in the console and/or in the log file.

        """

        formula = structure.composition.reduced_formula
        self.insert_data(
            functional+'.bulk_data', 
            {'mpid': mp_id,
             'formula': formula,
             'structure_fromMP': structure.as_dict()}, 
             message=message)
    
    def get_bulk_from_db(self, mp_id, functional, warning=False):
        """
        Get the data about bulk from the triboflow high level database.
        
        Parameters
        ----------
        mp_id : str
            Materials Project id of the structure.
        
        functional : str
            Functional. It could be PBE por SCAN.

        warning : bool
            Raise a warning instead of an error if the structure bulk is not
            found in the database. The default is False.
        
        Raises
        ------
        NavigatorError
            If there is no result from the query.

        Returns
        -------
        bulk : variable object
            Database object which contains the data for the selected bulk.

        """

        bulk = self.find_data(
            functional+'.bulk_data',
            {'mpid': mp_id})
        
        if bulk:
            return bulk
        else:
            message = 'No bulk material with MP-ID {} is found in the ' \
                      '{}.bulk_data collection.'.format(mp_id, functional)
            if warning:
                log.warning(message)
                return None
            else:
                raise StuctNavigatorError(message)

    def get_slab_from_db(self, mp_id, functional, miller, warning=False):
        """
        Get the data about slab from the triboflow high level database.

        Parameters
        ----------
        mp_id : str
            Materials Project id of the structure.
        
        functional : str
            Functional. It could be PBE por SCAN.
        
        miller : str
            Miller indices for the slab as a list of three integers.
        
        warning : bool
            Raise a warning instead of an error if the structure slab is not
            found in the database. The default is False.
        
        Raises
        ------
        StuctNavigatorError
            If there is no result from the query.

        Returns
        -------
        slab : variable object
            Database object which contains the data for the selected slab.
        
        """

        slab = self.find_data(
            functional+'.slab_data',
            {'mpid': mp_id, 'miller': miller})
        
        # Return the slab or alternatively raise an error or warning
        if slab:
            return slab
        else:
            message = 'No slab with MP-ID {} and miller indices ' \
                      '{} was found in the {}.slab_data collection.' \
                      .format(mp_id, miller, functional)
            if warning:
                log.warning(message)
                return None
            else:
                raise StuctNavigatorError(message)

    def get_interface_from_db(self, name, functional, warning=False):
        """
        Get the data about intreface from the eos database.

        Parameters
        ----------
        name : str
            Unique name of the interface as 
        
        functional : str
            Functional. It could be PBE por SCAN.

        warning : bool
            Raise a warning instead of an error if the selected interface is not
            found in the database. The default is False.

        Raises
        ------
        StuctNavigatorError
            If there is no result from the query.
        
        Returns
        -------
        interface : variable object
            Database object which contains the data for the selected slab.
    
        """
        
        interface = self.find_data(
            functional+'.interface_data',
            {'name': name})
        
        if interface:
            return interface
        else:
            message = 'No interface with name {} was found in the ' \
                      '{}.interface_data collection.'.format(name, functional)
            if warning:
                log.warning(message)
                return None
            else:
                raise StuctNavigatorError(message)

    def get_last_bmd_data_from_db(self, formula):
        """
        Get the last bulkmodule from the FireWorks database.

        A query is made to the eos collection in the Fireworks database for a
        given chemical formula. The results are sorted ascending by creation 
        date and than the last one is returned as a dictionary.

        Parameters
        ----------
        formula : str
            Chemical formula of the interface to query.

        Returns
        -------
        interface : variable object
            Interface generated by the get_wf_bulk_modulus workflow of atomate.

        """

        interface = self.db.eos.find(
            {'formula_pretty': formula}).sort('created_at', pymongo.DESCENDING)
        return interface[0]
    
class NavigatorMP:
    """
    This class is a high level interface for connecting with the Materials
    Project database (https://materialsproject.org/).

    These functions are taken from the TriboFlow package, Michael Wolloch.

    Attributes
    ----------
    mpr : database connection
        Connection to the Materials Project Database.

    Methods
    -------
    get_low_energy_structure(chem_formula, mp_id, print_info)
        Get the low energy structure from the Materials Project database.
    get_property_from_mp(mp_id, prop)
        Get the searched property from the Materials Project database.

    """

    def __init__(self):
        # Start a connection with the Materials Project database
        with MPRester() as mpr:
            self.__mpr = mpr

    def get_mpid_from_formula(self, chem_formula):
        """
        Return the MP_ID given the chemical formula.

        Parameters
        -----------
        chem_formula : str
            Chemical formula of the of the structure.
            e.g.: NaCl, Fe2O3, SiO, FeCW.
   
        Returns
        -------
        mp_id : str
            Materials Project ID of the structure.

        """
        mp_id = self.__mpr.query(criteria={'pretty_formula': chem_formula, 
                                      'e_above_hull': 0.0},
                                 properties=['material_id'])

        if len(mp_id) == 0 or mp_id is None:
            raise NavigatorMPError('{} has not been found in the Materials Project'
                                   ' database.'.format(chem_formula))
        return mp_id[0]['material_id']

    
    def get_low_energy_structure(self, chem_formula, mp_id=None, 
                                 print_info=False):
        """
        Search MaterialsProjects for structure.

        A function that searches the MaterialsProject Database
        for structures that match the given chemical formula
        and selects the one with the lowest formation energy
        per atom. If mp_id is given, the structure with that mp_id will
        be returned.

        Parameters
        ---------- 
        chem_formula : str
            Chemical formula of the structure.
            e.g.: NaCl, Fe2O3, SiO, FeCW.
        mp_id : str or None, optional
            Materials Project ID of the desired structure. The default is None.
            e.g.: 'mp-990448'.       
        print_info : bool or None, optional
            Whether to print some information about the collected structure.
            The default is False
             
        Returns
        -------
        struct : pymatgen.core.structure.Structure
            Tuple containing several information about the desired structure.

        mp_id : str
            Materials Project ID for the given chemical formula.)

        Examples
        --------
        Calling the method on 'NaCl' returns:
            (Structure Summary
                Lattice
                    abc : 4.024635423838785 4.024635423838785 4.024635423838785
                angles : 60.00000000000001 60.00000000000001 60.00000000000001
                volume : 46.09614833243692
                    A : 0.0 2.845847 2.845847
                    B : 2.845847 0.0 2.845847
                    C : 2.845847 2.845847 0.0
                PeriodicSite: Na (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
                PeriodicSite: Cl (2.8458, 2.8458, 2.8458) [0.5000, 0.5000, 0.5000])

        """

        if mp_id:
            struct = self.__mpr.get_structure_by_material_id(mp_id)

            return struct, mp_id
        else:
            id_list = self.__mpr.query(
                criteria={'pretty_formula': chem_formula,
                          'e_above_hull': 0.0},
                properties=['material_id'])

            if id_list == []:
                raise NavigatorMPError('{} has not been found in the MaterialsProject'
                                       'database'.format(chem_formula))
            else:
                mp_id = id_list[0]['material_id']
                struct = self.__mpr.get_structure_by_material_id(mp_id)

                return struct, mp_id

    def get_property_from_mp(self, mp_id, properties):
        """
        Get a certain property for a single material from the Materials Project 
        database.

        Parameters
        ----------
        mp_id : str
            Valid materials project ID.

        properties : list of str
            Properties for which to query. If just one property is not supported
            a warning will be issued.

        Returns
        -------
        searched_properties : dict
            The dictionary contains the all the searched properties.
            e.g.: mp_id='mp-990448' and properties=['energy', 'energy_per_atom'] 
                  the result is:
                  {'energy': -18.43845026, 'energy_per_atom': -9.21922513}

        """

        # Check if there is at least one not supported property
        supported_properties = self.__mpr.supported_task_properties
        not_found = [p in supported_properties for p in properties]
        if not all(not_found):
            print('In {} there is one or more not supported properties. \n'
                  'The supported properties are: \n {}.'
                  .format(properties, supported_properties))

            return None
        else:
            query = self.__mpr.query(criteria={'material_id': mp_id},
                                     properties=properties)
            if not query or len(query) == 0:
                print('The query has return nothing with the MP_ID: {}. \n'
                      'Please check if the MP_ID is a valid id or if there '
                      'is nothing for the searched property.'.format(mp_id))
            searched_property = query[0]

            return searched_property


# ============================================================================
# Functions to convert images to bytes and viceversa
# ============================================================================

def convert_bytes_to_image(bytes_object):
    """
    Convert an image saved as bytes for storing in MongoDB to an PIL image.

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

def convert_image_to_bytes(path_to_fig):
    """
    Convert an image to bytes for storage in MongoDB database.

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

def image_bytes_converter(data, to_image=True):
    """
    Convert an image to a bytes-object or a bytes-object to a PIL image. 
    Useful for storing data in MongoDB database.

    Parameters
    ----------
    data : bytes or str
        bytes-object or a path to an image.
        
    to_image : bool, optional
        Decide whether to convert to a PIL image or bytes. The default is True.

    Returns
    -------
    data_conv : PIL.PngImagePlugin.PngImageFile or bytes
        Converted data type

    """
    
    # We have bytes and convert to image
    if to_image:
        data_conv = convert_bytes_to_image(data)
      
    # We have (a path to an) image and convert to bytes
    else:
        data_conv = convert_image_to_bytes(data)
            
    return data_conv

