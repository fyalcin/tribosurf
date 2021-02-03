#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:02:48 2021

Classes to generate subworkflows for calculating slab properties.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

# ============================================================================
# Classes
# ============================================================================

class SlabWFs:
    """
    Collection of static methods to manipulate slab structures.
    
    """
    
    @staticmethod
    def converge_thickness_surfene():        
        """
        Method description
        """
        pass
        
    @staticmethod
    def converge_thickness_alat():        
        """
        Method description
        """
        pass

class SlabThicknessError(Exception): pass