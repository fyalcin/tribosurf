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

class GenerateSlabsErrors(Exception):
    """ Error in generating slabs.
    """
    @staticmethod
    def check_thickness(thickness):
        if not isinstance(thickness, (list, float, int)):
            raise GenerateSlabsErrors("Wrong type for argument: thickness. "
                                      "Allowed types: list, float, int.")
    
    @staticmethod
    def check_slabname(slab_name):
        if not isinstance(slab_name, (str, list)):
            raise GenerateSlabsErrors("Wrong type for argument: slab_name. "
                                      "Allowed types: list of str, str.")
        elif not all([isinstance(x, str) for x in list(slab_name)]):
            raise GenerateSlabsErrors("Wrong type for elements of list: "
                                       "slab_name. Allowed types: str.")            
            
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