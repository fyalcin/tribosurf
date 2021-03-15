#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:56:54 2021

Test the Firetasks of the `core` module.

Author: Omar Chehaimi
Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Omar Chehaimi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 22nd, 2021'

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.firetasks.slabs import FT_GenerateSlabs
from triboflow.firetasks.core import FT_RelaxStructure, FT_MoveTagResults
from triboflow.utils.database import Navigator, NavigatorMP

# Retriving the bulk structure from MP database and save it in our test database
chem_formula = 'Al'
mp_id = 'mp-134'

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
tag = 'relaxation_test'
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

# Move the relaxed structure to our custom database
ft_movetag = FT_MoveTagResults(
      mp_id=mp_id,
      collection_from=functional+collection,
      collection_to=functional+collection,
      db_file=db_file,
      database_from=None,
      database_to='test',
      miller=miller,
      tag=tag,
      check_entry=[
      ['relaxed_structure', 
      'data_'+chem_formula, 
      'calc_output']
      ],
      entry_to=[
      ['relaxed_structure', 
      'data_'+chem_formula, 
      'calc_output']
      ],
      entry_from=[
      ['output', 'structure'],
      ['nsites'],
      ['output', 'density'],
      ['output', 'energy'],
      ['output', 'energy_per_atom' ],
      ['output', 'bandgap'],
      ['output', 'forces'],
      ['output', 'stresses'],                                         
      ['_id']
      ],
      struct_kind='slab',
      override=True,
      cluster_params={})

wf = Workflow([Firework([ft_relax, ft_movetag])])

# Launch the calculation
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)