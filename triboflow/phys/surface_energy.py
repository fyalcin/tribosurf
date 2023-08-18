#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:00:54 2021

Functions to calculate the surface energy of a slab.

The module contains the following functions:
    - calculate_surface_energy

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Gabriele Losi"
__copyright__ = (
    "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
)
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 22nd, 2021"

import numpy as np
from pymatgen.core.structure import Structure


def calculate_surface_energy(output_list, sym_surface=True):
    """
    Calculate the surface energy by passing a list containing all the dictionary
    with the output data. The first element is treated to be the bulk oriented
    along a specific miller index direction.

    Parameters
    ----------
    output_list : list of dict
        List of dictionaries where each one contains the output of a VASP
        simulation. Basically output_list can be whatever dictionary, fundamental
        keys that should be present are: structure, energy, energy_per_atom.

    sym_surface : bool, optional
        If the surfaces have to be considered symmetric or not.

    Returns
    -------
    surfene : list of floats
        List containing the surface energies calculated from output_list.
    """

    # Take out the energy per atom of the bulk
    energy_per_atom = output_list[0]["energy_per_atom"]

    # Calculate the surface area of a slab
    bulk_latvecs = Structure.from_dict(output_list[1]["structure"]).lattice.matrix
    area = np.linalg.norm(np.cross(bulk_latvecs[0], bulk_latvecs[1]))

    # Loop over the slab elements of output_list and calculate surface energy
    surfene = np.array([])
    for d in output_list[1:]:
        energy = d["energy"]
        nsites = d["nsites"]
        surfene = np.append(
            surfene, 16.02176565 * (energy - energy_per_atom * nsites) / area
        )

    # Divide the surface energies by two if the surfaces are symmetric
    if sym_surface:
        surfene /= 2.0

    return np.array(surfene)
