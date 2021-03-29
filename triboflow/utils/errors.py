#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:03:22 2021

Custom errors for any Firetasks and Workflow.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 8th, 2021'


class GeneralErrorFT(Exception):

    @staticmethod
    def check_collection(collection, functional = ['PBE', 'SCAN'], 
                         data = ['bulk_data', 'slab_data', 'interface_data']):
        
        allowed_collections = [f + '.' + d for f in functional for d in data]

        if collection not in allowed_collections:
            raise RelaxStructureError('Wrong value for collection. Allowed '
                                        'values: "bulk_data", "slab_data", '
                                        '"interface_data"')

class SlabOptThickError(GeneralErrorFT): 
    """ Error in Slab Optimal Thickness workflows and dependencies.
    """
    pass

class GenerateSlabsError(GeneralErrorFT):
    """ Error in generating slabs.
    """

    @staticmethod
    def check_miller(miller):
        if not isinstance(miller, list):
            raise GenerateSlabsError("Wrong type for argument: miller. "
                                     "Allowed types: list")
        elif any([isinstance(x, list) for x in miller]) and not all([isinstance(x, list) for x in miller]):
            raise GenerateSlabsError("Wrong type for elements of list: "
                                     "miller. Allowed types: list.")

    @staticmethod
    def check_thickness(thickness):
        if not isinstance(thickness, (list, float, int)):
            raise GenerateSlabsError("Wrong type for argument: thickness. "
                                     "Allowed types: list, float, int.")
        if isinstance(thickness, list):
            if not all([isinstance(x, (float, int)) for x in thickness]):
                raise GenerateSlabsError("Wrong type for elements of list: "
                                         "miller. Allowed types: int, float.")

    @staticmethod
    def check_vacuum(vacuum):
        if not isinstance(vacuum, (list, float, int)):
            raise GenerateSlabsError("Wrong type for argument: vacuum. "
                                     "Allowed type: list")
        if isinstance(vacuum, list):
            if not all([isinstance(x, (float, int)) for x in vacuum]):
                raise GenerateSlabsError("Wrong type for elements of list: "
                                         "vacuum. Allowed types: list.")

    @staticmethod
    def check_entry(entry):
        if not isinstance(entry, (str, list)):
            raise GenerateSlabsError("Wrong type for argument: entry. "
                                     "Allowed types: list of str, str.")

    @staticmethod
    def check_inputs(miller, thickness, vacuum, entry):
        
        GenerateSlabsError.check_miller(miller)
        GenerateSlabsError.check_thickness(thickness)
        GenerateSlabsError.check_vacuum(vacuum)
        GenerateSlabsError.check_entry(entry)

        # If one is a list, both of them should be.
        if isinstance(thickness, list):
            if not isinstance(entry, list):
                raise GenerateSlabsError("Wrong type for arguments: thickness, "
                                         "entry. If one is a list, both "
                                         "should be.")
            if len(thickness) != len(entry):
                raise GenerateSlabsError("Wrong arguments: thickness, entry. " 
                                         "They should have the same length if lists")               

class RelaxStructureError(GeneralErrorFT):
    """ Error when running a relax calculation
    """

    @staticmethod
    def is_data(structure, mp_id, functional):
        if structure is None:
            raise RelaxStructureError('No entry found in DB for a '
                                      'structure with mpid: {}, functional: {}'
                                      .format(mp_id, functional))

class MoveTagResultsError(GeneralErrorFT):
    """ Errors when moving data between two different db_file/database/collection.
    """
    
    @staticmethod
    def check_entry(entry, msg='entry'):

        if isinstance(entry, list):
            if not all([isinstance(n, list) for n in entry]):
                raise MoveTagResultsError("Wrong type for arguments: {}. "
                                          "Allowed type: list".format(msg))
        
        # elif not isinstance(entry, str):
        #     raise MoveTagResultsError("Wrong argument: {}. "
        #                               "Allowed types: str, list".format(msg))

    @staticmethod
    def check_entries(entry_check, entry_to, entry_from):

        MoveTagResultsError.check_entry(entry_check, 'entry_check')
        MoveTagResultsError.check_entry(entry_to, 'entry')
        MoveTagResultsError.check_entry(entry_from, 'entry_tag')

        if len(list(entry_to)) != len(list(entry_from)):
            raise MoveTagResultsError("Wrong arguments: entry or entry_tag. If "
                                      "converted to str, should have same length")

class SurfaceEnergyError(GeneralErrorFT):
    """ Error in surface energy Firetask.
    """
    pass

class SubWFError(Exception):
    """ Error in reading the subworkflow parameters.
    """
    pass

class NavigatorError(Exception):
    """ Errors in database initialization and queries, for Navigator class.
    """
    pass

class StuctNavigatorError(NavigatorError):
    """ Errors in queries for structures handling in low and high level DBs.
    """ 
    pass

class NavigatorMPError(Exception):
    """ Errors in querying the online MP database.
    """
    pass

class ReadParamsError(Exception):
    """ Errors when reading, parsing or retrieving data or parameters in FTs.
    """
    pass

class WriteParamsError(Exception):
    """ Errors when writing data to dictionary to be stored in DB.
    """
    pass