#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:27:12 2021

Test the Firetasks of the `slabs` module.

Author: Gabriele Losi (glosi000)
Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 22nd, 2021'

from pymatgen.core.structure import Structure
from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.firetasks.slabs import FT_GenerateSlabs
from triboflow.utils.database import StructureNavigator


# Define input parameters
mp_id = 'mp-126'
functional = 'PBE'
miller = [[1, 0, 0]] * 7
db_file = None
high_level = 'tribchem'
database = high_level
thickness = [0, 4, 6, 8, 10, 12, 14]
thick_bulk = 4
vacuum = 10
ext_index = 0
in_unit_planes = True
entry = [['thickness', 'data_' + str(thk), 'unrelaxed'] for thk in thickness]

# Call the Navigator and retrieve the bulk structure
nav = StructureNavigator(db_file, high_level)
structure = nav.get_bulk_from_db(mp_id, functional, warning=True)['structure_fromMP']
structure = Structure.from_dict(structure)

# Instantiate the Firetask and create a WF
ft = FT_GenerateSlabs(structure=structure,
                      mp_id=mp_id,
                      miller=miller,
                      collection=functional+'.slab_data',
                      db_file=db_file,
                      database=database,
                      thickness=thickness,
                      thick_bulk=thick_bulk,
                      vacuum=vacuum,
                      symmetrize=False,
                      ext_index=ext_index,
                      in_unit_planes=in_unit_planes,
                      entry=entry)
wf = Workflow([Firework([ft])])

# Run the workflow
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)
