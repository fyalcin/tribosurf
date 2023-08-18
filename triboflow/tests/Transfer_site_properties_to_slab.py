#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 13:01:38 2020

@author: mwo
"""

import numpy as np
from pymatgen import MPRester
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from pymatgen.analysis.structure_matcher import *

NiO_struct = MPRester().get_structure_by_material_id("mp-19009")
Fe_struct = MPRester().get_structure_by_material_id("mp-13")


def make_magnetic_standard_structure(struct):
    conv_struct = sga(struct).get_conventional_standard_structure(keep_site_properties=True)
    my_lattice = struct.lattice
    conv_lattice = conv_struct.lattice

    # Get the tranisition matrix from the initial to the conventional lattice
    TM_MyToConv = np.dot(my_lattice.inv_matrix, conv_lattice.matrix)

    # This scaling factor is not working generally! It will approach inf if one
    # or more entries in the TM are close to 0!
    scale_factor = 1 / np.amin(np.absolute(TM_MyToConv))
    SC_matrix = TM_MyToConv * scale_factor

    my_supercell = struct.copy()
    my_supercell.make_supercell(SC_matrix)

    conv_supercell = conv_struct.copy()
    conv_supercell.make_supercell(scale_factor)

    matcher = StructureMatcher(primitive_cell=False)
    matched_list = matcher.get_mapping(my_supercell, conv_supercell)

    for i, site in enumerate(conv_supercell.sites):
        mag_moment = my_supercell.site_properties["magmom"][matched_list[i]]
        conv_supercell.replace(
            i, conv_supercell.species[i], properties={"magmom": mag_moment}
        )
    return conv_supercell


# Fe_conv = make_magnetic_standard_structure(Fe_struct)
NiO_conv = make_magnetic_standard_structure(NiO_struct)
