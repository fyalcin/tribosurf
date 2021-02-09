#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:03:22 2021

Custom errors for any Firetasks and Workflow.

@author: glosi000
"""

class SlabOptThickError(Exception): 
    """ Error in Slab Optimal Thickness workflows and dependencies
    """
    pass

class GenerateSlabsErrors(Exception):
    """ Error in generating slabs
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
        elif not all([isinstance(x, str) for x in list(slab_name]):
            raise GenerateSlabsErrors("Wrong type for elements of list: "
                                       "slab_name. Allowed types: str.")            



                

class ReadSubWFsError(Exception):
    """ Error in reading the subworkflow parameters
    """
    pass
