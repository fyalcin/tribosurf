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

class ReadSubWFsError(Exception):
    """ Error in reading the subworkflow parameters
    """
    pass
