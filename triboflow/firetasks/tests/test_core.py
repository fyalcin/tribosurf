#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:56:54 2021

Test the Firetasks of the `core` module.

Author: Gabriele Losi (glosi000)
Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.firetasks.core import FT_MoveTagResults
from triboflow.utils.database import Navigator


# Define input parameters
mp_id = 'mp-126'
collection_from = "tasks"
collection_to = "PBE.slab_data"
tag = True
tag_key = "transfer_test"
db_file = None
database_from = None
database_to = 'test'
miller = [0, 1, 2]
check_entry = None
entry_to = [['test', 'energy'], ['test', 'energy2']]
entry_from = [['data', 'energy'], ['data', 'energy2']]
override = False
cluster_params = {}

# Initialize the collections for the tests
nav = Navigator(db_file=db_file)
data = {'test': 
           {'energy1': 1,
            'energy2': 2,
            'energy3': 3,
            }
}
nav.insert_data(collection='tasks', data={'energy': 100})

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
                       check_entry=check_entry,
                       entry_to=entry_to,
                       entry_from=entry_from,
                       override=override,
                       cluster_params=cluster_params)
                       
wf = Workflow([Firework([ft])])

# Run the workflow
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)
