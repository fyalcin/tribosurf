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

import os
from pathlib import Path, PurePosixPath

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.core.structure import Structure
from triboflow.utils.database import NavigatorMP, StructureNavigator
from triboflow.workflows.slabs_wfs import SlabWF

# Get the bulk
# formula = 'Cu'
# mid = 'mp-30'
# nav_mp = NavigatorMP()
formula = "Al"
mid = "mp-134"
nav = StructureNavigator("auto", "triboflow")
bulk_dict = nav.get_bulk_from_db(mid, "PBE")
structure = Structure.from_dict(bulk_dict["structure_equiVol"])
comp_params = bulk_dict["comp_parameters"]


# structure, mid = nav_mp.get_low_energy_structure(
#     chem_formula=formula,
#     mp_id=mid)

# Surface generation tests
wf = SlabWF.conv_slabthick_surfene(
    structure=structure,
    mp_id=mid,
    comp_params=comp_params,
    miller=[1, 1, 1],
    thick_min=3,
    thick_max=18,
    thick_incr=1,
    parallelization="low",
    vacuum=25,
    conv_thr=0.01,
)

# Launch the calculation
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
# rapidfire(lpad)
