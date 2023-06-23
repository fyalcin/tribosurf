#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 4 11:36:15 2021

Tests for surface energy calculation.
The following tests are about the slabs generation and the following surface 
energy calculation for an aluminum system.

Author: Omar Chehaimi
Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Omar Chehaimi"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "March 4th, 2021"


from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from triboflow.utils.database import Navigator, NavigatorMP, StructureNavigator
from triboflow.workflows.surfene_wfs import SurfEneWF


# Get the bulk from the online Database: Materials Project
nav = StructureNavigator("auto", True)
mpid = "mp-134"
mat_dict = nav.get_bulk_from_db(functional="PBE", mp_id=mpid)
structure = Structure.from_dict(mat_dict["primitive_structure"])
# formula = 'Mg'
# mid = 'mp-110'
# nav_mp = NavigatorMP()
# structure, mid = nav_mp.get_low_energy_structure(
#    chem_formula=formula,
#    mp_id=mid)

# Get the bulk from a local simple Poscar
# structure = Structure.from_file('POSCAR')
# mid = 'custom-1'

# Surface generation tests
wf = SurfEneWF.conv_surface_energy(
    structure=structure,
    mp_id=mpid,
    miller=[0, 0, 1],
    thick_min=2,
    thick_max=8,
    thick_incr=2,
    parallelization=None,
    add_static=True,
    db_file="auto",
    low_level="FireWorks",
    high_level=True,
)

# Launch the calculation
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
# rapidfire(lpad)

# Check if the structures are correct
# a = Navigator()
# b = a.find_data('PBE.slab_data', {'mpid': 'mp-110'})
# oriented_bulk = Structure.from_dict(b['thickness']['data_0']['input'])
# min_thick_slab = Slab.from_dict(b['thickness']['data_2']['input'])
# max_thick_slab = Slab.from_dict(b['thickness']['data_3']['input'])
