#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:00:47 2021

Functions to manage crystalline structures and solid state physics.

The module contains the following functions:
    - orient_bulk
    - generate_slabs

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Gabriele Losi"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 22nd, 2021"

from pymatgen.core.surface import SlabGenerator


# ============================================================================
# Functions to deal with crystalline slabs
# ============================================================================


def orient_bulk(
    structure,
    miller,
    thickness,
    primitive=False,
    lll_reduce=False,
    in_unit_planes=True,
):
    """
    Orient a bulk unit cell along a direction identified by Miller indexes.

    """

    # Generate the oriented bulk
    slabgen = SlabGenerator(
        initial_structure=structure,
        miller_index=miller,
        primitive=primitive,
        lll_reduce=lll_reduce,
        in_unit_planes=in_unit_planes,
        min_slab_size=thickness,
        min_vacuum_size=0,
    )

    bulk_miller = slabgen.oriented_unit_cell

    return bulk_miller
