#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 17:58:50 2021

Test the generation of the computational parameters to be used by VASP.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Gabriele Losi"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "7th April, 2021"


from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from triboflow.utils.tests.vasp_tools_debug import (
    get_custom_vasp_relax_settings,
)


bulk = Structure.from_file("../../workflows/tests/POSCAR")

slabgen = SlabGenerator(
    initial_structure=bulk,
    miller_index=[0, 0, 1],
    center_slab=True,
    primitive=False,
    lll_reduce=True,
    in_unit_planes=True,
    min_slab_size=2,
    min_vacuum_size=10,
)

# Select the ext_index-th slab from the list of possible slabs
s = slabgen.get_slabs()
# bonds=Nol,
#                       ftol=ftol,
#                       tol=tol,
#                       repair=repair,
#                       max_broken_bonds=max_broken_bonds,
#                       symmetrize=symmetrize)
s = s[0]

comp_parameters = {}

vis_1 = get_custom_vasp_relax_settings(
    bulk, comp_parameters, "bulk_shape_relax"
)
vis_2 = get_custom_vasp_relax_settings(s, comp_parameters, "slab_shape_relax")
