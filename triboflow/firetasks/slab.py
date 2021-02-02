#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:09:03 2021

Firetasks to converge the slab thickness either calculating the surface 
energy or the lattice parameter for structures with increasing number of 
atomic layers.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase
from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator, NavigatorMP

# ============================================================================
# Firetasks
# ============================================================================

@explicit_serialize
class FT_SlabThicknessConvo(FiretaskBase):
    """
    Firetask description...
    
    """
    
    _fw_name = 'Slab Thickness convergence'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file']
    
    def run_task(self, fw_spec):
        pass


# ============================================================================
# Functions
# ============================================================================