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


class SlabOptThickError(Exception): 
    """ Error in Slab Optimal Thickness workflows and dependencies.
    """
    pass

class GenerateSlabsError(Exception):
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
            if not all([isinstance(x, float, int) for x in thickness]):
                raise GenerateSlabsError("Wrong type for elements of list: "
                                         "miller. Allowed types: int, float.")

    @staticmethod
    def check_vacuum(vacuum):
        if not isinstance(vacuum, (list, float, int)):
            raise GenerateSlabsError("Wrong type for argument: vacuum. "
                                     "Allowed type: list")
        if isinstance(check_vacuum, list):
            if not all([isinstance(x, (float, int)) for x in vacuum]):
                raise GenerateSlabsError("Wrong type for elements of list: "
                                         "vacuum. Allowed types: list.")

    @staticmethod
    def check_slabname(slab_name):
        if not isinstance(slab_name, (str, list)):
            raise GenerateSlabsError("Wrong type for argument: slab_name. "
                                     "Allowed types: list of str, str.")

    @staticmethod
    def check_inputs(miller, thickness, vacuum, slab_name):
        
        GenerateSlabsError.check_miller(miller)
        GenerateSlabsError.check_thickness(thickness)
        GenerateSlabsError.check_vacuum(vacuum)
        GenerateSlabsError.check_slab_name(slab_name)

        # If one is a list, both of them should be.
        if isinstance(thickness, list):
            if not isinstance(slab_name, list):
                raise GenerateSlabsError("Wrong type for arguments: thickness, "
                                         "slab_name. If one is a list, both "
                                         "should be.")
            if len(thickness) != len(slab_name):
                raise GenerateSlabsError("Wrong arguments: thickness, slab_name. " 
                                         "They should have the same length if lists")               

class RelaxStructureError(Exception):
    """ Error when running a relax calculation
    """

    @staticmethod
    def check_collection(collection):
        if collection not in ['bulk_data', 'slab_data', 'interface_data']:
            raise RelaxStructureError('Wrong value for collection. Allowed '
                                      'values: "bulk_data", "slab_data", '
                                      '"interface_data"')

    @staticmethod
    def is_data(structure, mp_id, functional):
        if structure is None:
            formula = structure.composition.reduced_formula
            raise RelaxStructureError('No entry found in DB {} for a '
                                      'structure with mpid: {}, functional: {}'
                                      .format(formula, mp_id, functional))

class MoveStructInDBError(RelaxStructureError):
    """ Errors when moving data between two different db_file/database/collection.
    """
    
    @staticmethod
    def check_name(name, msg='name'):

        if isinstance(name, list):
            if not all([isinstance(n, list) for n in name_tag]):
                raise MoveStructInDBError("Wrong type for arguments: {}. "
                                          "Allowed type: list".format(msg))
        
        elif not isinstance(name, str):
            raise MoveStructInDBError("Wrong argument: {}. "
                                      "Allowed types: str, list".format(msg))

    @staticmethod
    def check_name_str(name, msg='name'):
        if not isinstance(name, str):
            raise MoveStructInDBError("Wrong argument: {}. Allowed types: str"
                                      "".format(msg))

    @staticmethod
    def check_names(name_check, name, name_tag):

        MoveStructInDBError.check_name_str(name_check, 'name_check')
        MoveStructInDBError.check_name(name, 'name')
        MoveStructInDBError.check_name(name_tag, 'name_tag')

        if len(list(name)) != len(list(name_tag)):
            raise MoveStructInDBError("Wrong arguments: name or name_tag. If "
                                      "converted to str, should have same lenght")

class SurfaceEnergyError(Exception):
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

class ReadWriteParamsError(Exception):
    """ Errors when writing data to dictionary to be stored in DB.
    """
    pass