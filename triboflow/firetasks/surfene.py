#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:54:53 2021

Firetasks to calculate the surface energy for a bulk and a given orientation.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase
from atomate.utils.utils import env_chk
from pymatgen.transformations.advanced_transformations import SlabTransformation

from triboflow.utils.database import Navigator, NavigatorMP


# ============================================================================
# Firetasks
# ============================================================================

@explicit_serialize
class FT_SurfaceEnergy(FiretaskBase):
    """
    Firetask description...
    
    """
    
    _fw_name = 'Surface Energy calculation'
    required_params = ['bulk', 'slab', 'tag']
    optional_params = ['db_file']
    
    def run_task(self, fw_spec):
        



# ============================================================================
# Functions
# ============================================================================
