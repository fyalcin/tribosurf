#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:27:12 2021

The the Firetasks of the `slabs` module.

Author: Gabriele Losi (glosi000)
Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna
"""

from triboflow.firetasks.slabs import FT_GenerateSlabs
from triboflow.utils.database import StructureNavigator
import os

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from pymatgen.core.surface import Structure

from fireworks import Firework, Workflow

currentdir = os.path.dirname(__file__)

mp_id = 'mp-126'
functional = 'PBE'
miller = [[1, 0, 0], [3, 3, 3], [1, 2, 2], [11, 0, 0]]
db_file = None
high_level = "tribchem"
database = high_level
thickness = [12,12,12,12]
thick_max = 14
vacuum=[10, 10, 10, 10]
ext_index=0
in_unit_planes=True
slab_name=['unrelaxed', 'unrelaxed', 'unrelaxed', 'unrelaxed']

nav = StructureNavigator(db_file, high_level)
structure = nav.get_bulk_from_db(mp_id, functional, warning=True)['structure_fromMP']
structure = Structure.from_dict(structure)

ft = FT_GenerateSlabs(structure=structure,
                      mp_id=mp_id,
                      miller=miller,
                      collection=functional+'.slab_data',
                      db_file=db_file,
                      database=database,
                      thickness=thickness,
                      thick_bulk=thick_max,
                      vacuum=vacuum,
                      symmetrize=False,
                      ext_index=ext_index,
                      in_unit_planes=in_unit_planes,
                      entry=slab_name)

wf = Workflow([Firework([ft])])

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)