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
from fireworks import Workflow, Firework

from triboflow.utils.database import NavigatorMP
from triboflow.workflows.slabs_wfs import SlabWF
from triboflow.firetasks.run_slabs_wfs import FT_SlabOptThick


# Get the bulk from the online Database: Materials Project
formula = 'Cu'
functional = 'PBE'
miller=[0, 0, 1]
mid = 'mp-30'
nav_mp = NavigatorMP()
structure, mid = nav_mp.get_low_energy_structure(
   chem_formula='formula', 
   mp_id=mid)

# Get the bulk from a local simple Poscar
# structure = Structure.from_file('POSCAR')
# mid = 'custom-1'


# FIRST TEST, RUN THE SLAB OPT THICKNESS AS A SUBWORKFLOW
# ft = FT_SlabOptThick(mp_id=mid,
#                      miller=miller, 
#                      functional=functional,
#                      thick_min=3, 
#                      thick_max=6,
#                      thick_incr=1,
#                      db_file=None,
#                      low_level=None,
#                      high_level='triboflow',
#                      vacuum=10,
#                      in_unit_planes=True,
#                      ext_index=0,
#                      conv_thr=0.01,
#                      parallelization='low',
#                      cluster_params={},
#                      override=True)
# fw = Firework(ft, name = 'Start a subworkflow with a single FT')
# wf_1 = Workflow([fw], name = 'Start a subworkflow to converge slab thickness')


# SECOND TEST, RUN THE OPTIMIZATION WORKFLOWS ALONE
wf_2 = SlabWF.conv_slabthick_surfene(structure=structure,
                                      mp_id=mid,
                                      miller=miller,
                                      thick_min=3, 
                                      thick_max=6,
                                      thick_incr=1,
                                      db_file=None,
                                      low_level=None,
                                      high_level='triboflow',
                                      vacuum=10,
                                      in_unit_planes=True,
                                      ext_index=0,
                                      conv_thr=0.00001,
                                      parallelization='high',
                                      recursion=0,
                                      cluster_params={},
                                      override=False)

lpad = LaunchPad.auto_load()

# Run the second test
lpad.add_wf(wf_2)

# Run first test
# lpad.add_wf(wf_1)

rapidfire(lpad)
