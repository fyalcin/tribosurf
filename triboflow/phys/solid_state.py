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
__copyright__ = (
    "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
)
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 22nd, 2021"

from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from hitmen_utils.shaper import Shaper
from triboflow.utils.structure_manipulation import transfer_average_magmoms


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


def generate_slabs(
    structure,
    miller,
    thickness,
    vacuum,
    thick_bulk=12,
    center_slab=True,
    primitive=True,
    lll_reduce=True,
    in_unit_planes=True,
    ext_index=0,
    bonds=None,
    ftol=0.1,
    tol=0.1,
    repair=False,
    max_broken_bonds=0,
    symmetrize=False,
):
    """
    Create and return a single slab or a list of slabs out of a structure.
    To return a list of slabs you need to provide `miller` as a list of lists,
    while `vacuum` and `thickness` can be either a single value (if you want all
    the slabs to have same thickness and vacuum) or a list of values.
    If `thickness` is 0, an oriented bulk is returned. The other values are
    passed as arguments to `SlabGenerator` to build the slabs.

    Examples
    --------
        miller : [0, 0, 1] or [[0, 0, 1], [1, 1, 1], [1, 1, 0]].

    """

    # Manage the arguments type in order to have lists
    if isinstance(miller, list) and not all([isinstance(m, list) for m in miller]):
        miller = [miller]
    if not isinstance(thickness, list):
        thickness = [thickness]
    if not isinstance(vacuum, list):
        vacuum = [vacuum]

    # Manage the length of the lists
    n = len(thickness)
    if len(vacuum) != n:
        vacuum *= n
    if len(miller) != n:
        miller *= n
    # SlabGenerator expects conventional unit cell, so we convert the structure accordingly.
    # As a result, we require input structure to be the primitive standard structure.
    conv_structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure(keep_site_properties=True)
    structure = transfer_average_magmoms(structure, conv_structure)

    slabs = []
    for hkl, thk, vac in zip(miller, thickness, vacuum):
        # If thk is zero, then we want to construct an oriented unit bulk by
        # our conversion, so we define a "fake" thickness > 1 to be used to
        # build the slabs, otherwise errors are raised by pymatgen
        thk_gen = 1 if thk == 0 else thk

        slabgen = SlabGenerator(
            initial_structure=structure,
            miller_index=hkl,
            center_slab=center_slab,
            primitive=primitive,
            lll_reduce=lll_reduce,
            max_normal_search=max([abs(index) for index in hkl]),
            in_unit_planes=in_unit_planes,
            min_slab_size=thk_gen,
            min_vacuum_size=vac,
        )

        s = slabgen.get_slabs(
            bonds=bonds,
            ftol=ftol,
            tol=tol,
            repair=repair,
            max_broken_bonds=max_broken_bonds,
            symmetrize=symmetrize,
        )

        # Case of an oriented bulk
        if thk == 0:
            s = s[ext_index].oriented_unit_cell.get_primitive_structure(
                constrain_latt={
                    "a": s[ext_index].lattice.a,
                    "b": s[ext_index].lattice.b,
                    "gamma": s[ext_index].lattice.gamma,
                }
            )

        # Case of a slab
        else:
            ouc = s[ext_index].oriented_unit_cell
            max_layer_spacing = max(Shaper.get_layer_spacings(ouc, ftol))
            s = [
                Shaper.resize(
                    struct=slab,
                    slab_thickness=thk,
                    vacuum_thickness=vac,
                    min_vac=1.2 * max_layer_spacing,
                )
                for slab in s
            ]
            s = s[ext_index]

        slabs.append(s)

    return slabs
