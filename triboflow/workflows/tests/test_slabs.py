#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 16:41:17 2021

Test the worflow to calculate the optimal slab thickness by surface energy.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'April 6th, 2021'


from pymatgen.core.structure import Structure
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.workflows.slabs_wfs import SlabWF
from triboflow.utils.database import NavigatorMP


# Get the bulk from the online Database: Materials Project
#formula = 'Cu'
#mid = 'mp-30'
#nav_mp = NavigatorMP()
#structure, mid = nav_mp.get_low_energy_structure(
#    chem_formula=formula, 
#    mp_id=mid)

# Get the bulk from a local simple Poscar
structure = Structure.from_file('POSCAR')
mid = 'custom-1'

wf = SlabWF.conv_slabthick_surfene(structure=structure,
                                   mp_id=mid,
                                   miller=[0, 0, 1],
                                   thick_min=2, 
                                   thick_max=10,
                                   thick_incr=2,
                                   db_file=None,
                                   low_level=None,
                                   high_level='triboflow',
                                   vacuum=10,
                                   in_unit_planes=True,
                                   ext_index=0,
                                   conv_thr=0.01,
                                   parallelization='low',
                                   recursion=False,
                                   cluster_params={},
                                   override=False)

# Launch the calculation
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
# rapidfire(lpad)
