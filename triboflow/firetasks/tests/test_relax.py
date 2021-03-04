#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:56:54 2021

Test the Firetasks of the `core` module.

Author: Omar Chehaimi
Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.firetasks.slabs import FT_GenerateSlabs
from triboflow.firetasks.core import FT_RelaxStructure
from triboflow.utils.database import Navigator, NavigatorMP

# Retriving the bulk structure from MP database and save it in our test database
chem_formula='Si'
mp_id = 'mp-149'

nav_mp = NavigatorMP()
structure, mp_id = nav_mp.get_low_energy_structure(chem_formula=chem_formula, 
                                                   mp_id=mp_id)

db_file = None
nav = Navigator(db_file=db_file, high_level='test')
nav.insert_data(collection='PBE.bulk_data', 
                data={'mpid': mp_id,
                      'formula': chem_formula,
                      'structure_fromMP': structure.as_dict(),
                      'comp_params': {}})

# Define input parameters for relaxing a bulk
functional = 'PBE'
collection = 'PBE.bulk_data'
entry = 'structure_fromMP'
tag = 'relaxation-test'
database = 'test'
relax_type = 'slab_pos_relax'
comp_params = {}
miller = [0, 1, 2]
check_key = None

# Call the relaxation
ft_relax = FT_RelaxStructure(mp_id=mp_id,
                             functional=functional,
                             collection=collection,
                             entry=entry,
                             tag=tag,
                             database=database,
                             relax_type=relax_type)

wf = Workflow([Firework([ft_relax])])

# Launch the calculation
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)