#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:56:54 2021

The the Firetasks of the `core` module.

Author: Gabriele Losi (glosi000)
Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.firetasks.core import FT_MoveTagResults


# Define input parameters
mp_id = 'mp-126'
collection_from = "tasks"
collection_to = "PBE.slab_data"
tag = True
tag_key = "transfer_test"
db_file = None
database_from = None
database_to = 'tribchem'
miller = [0, 1, 2]
entry_check = None
entry_to = [['test', 'energy'], ['test', 'energy2'], ['test', 'energy3']]
entry_from = [['data', 'energy'], ['data', 'energy2'], ['data', 'energy']]
override = False
cluster_params = {}

# Instantiate the Firetask and create a WF
ft = FT_MoveTagResults(mp_id=mp_id,
                       collection_from=collection_from,
                       collection_to=collection_to,
                       tag=tag,
                       tag_key=tag_key,
                       db_file=db_file,
                       database_from=database_from,
                       database_to=database_to,
                       miller=miller,
                       entry_check=entry_check,
                       entry_to=entry_to,
                       entry_from=entry_from,
                       override=override,
                       cluster_params=cluster_params)
                       
wf = Workflow([Firework([ft])])

# Run the workflow
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)
