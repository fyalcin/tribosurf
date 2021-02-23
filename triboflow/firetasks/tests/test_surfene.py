#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:58:05 2021

Test the Firetasks of the `surfene` module.

Author: Gabriele Losi (glosi000)
Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.firetasks.surfene import FT_SurfaceEnergy


# Define input parameters
mp_id = 'mp-126'
collection = "PBE.slab_data"
miller = [1, 0, 0]
thickness = [0, 4, 6, 8, 10, 12, 14]
entry = [['thickness', 'data_' + str(thk)] for thk in thickness]
db_file = None
database = "tribchem"

# Instantiate the Firetask and create a WF
ft = FT_SurfaceEnergy(mp_id=mp_id,
                      collection=collection,
                      miller=miller,
                      entry=entry,
                      db_file=db_file,
                      database=database)

wf = Workflow([Firework([ft])])

# Run the workflow
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)
