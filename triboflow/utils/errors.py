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
    def check_thickness(thickness):
        if not isinstance(thickness, (list, float, int)):
            raise GenerateSlabsError("Wrong type for argument: thickness. "
                                     "Allowed types: list, float, int.")
    
    @staticmethod
    def check_slabname(slab_name):
        if not isinstance(slab_name, (str, list)):
            raise GenerateSlabsError("Wrong type for argument: slab_name. "
                                     "Allowed types: list of str, str.")
        elif not all([isinstance(x, str) for x in list(slab_name)]):
            raise GenerateSlabsError("Wrong type for elements of list: "
                                     "slab_name. Allowed types: str.")

    @staticmethod
    def check_slabname_thickness(thickness, slab_name):

        GenerateSlabsError.check_thickness(thickness)
        GenerateSlabsError.check_slab_name(slab_name)

        # If one is a list, both of them should be.
        if isinstance(thickness, list):
            if not isinstance(slab_name, list):
                raise GenerateSlabsError("Wrong type for arguments: thickness, "
                                         "slab_name. They should be the")
            
class ReadSubWFsError(Exception):
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
