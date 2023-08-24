#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 11:57:19 2021

This module contains Firetasks useful to load atomic structures from local source
or online database, providing a starting point for any more complex workflow
initialized from scratch.

The module contains:

** Firetasks **:
    - FTCheckInput
    Checks a dictionary for essential keys and adds default values.

** Functions **:
    - unbundle_input
    - material_from_mp
    - read_input_dict
    - dict_consistency

    Author: Gabriele Losi (glosi000)
    Credits: The code is partially based on the original work of Michael 
    Wolloch, Triboflow package, University of Vienna
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Gabriele Losi"
__credits__ = "This module is based on the work of Michael Wolloch, TriboFlow"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "January 13th, 2021"

import os

from fireworks import FWAction, FiretaskBase, explicit_serialize

from triboflow.utils.mp_connection import MPConnection
from triboflow.utils.utils import load_defaults

currentdir = os.path.dirname(__file__)


# ============================================================================
# FireTasks
# ============================================================================


@explicit_serialize
class FTCheckInput(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values.

    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.

    Parameters
    ----------
    input_dict : list of str
        Location of the input parameters dictionary to be checked in the spec.

    read_key : str
        Dictionary key to be used to extract the corresponding list of keys
        and default values to compare the input to.

    output_dict_name: list of str, optional
        Location of the output dictionary that is going to be put into the
        spec. The default is 'out_params'.

    Returns
    -------
    dict
        A dictionary is pushed to the spec of the workflow. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys. The name in the spec is output_dict_name.

    """

    _fw_name = "Check input dictionary"
    required_params = ["input_dict", "read_key"]
    optional_params = ["output_dict_name", "fw_name"]

    def run_task(self, fw_spec):
        """Run the FireTask."""

        # Read the required and optional input parameters of the Firetask
        input_dict = self["input_dict"]
        output_dict_name = self.get("output_dict_name", "out_params")
        read_key = self["read_key"]

        # Create the output dictionary to be stored in DB
        out_dict = read_input_dict(input_dict, read_key)

        # mpid of minimum energy structure is used as default for materials
        if "mpid" in out_dict.keys() and out_dict["mpid"] is None:
            mp_connection = MPConnection()
            out_dict["mpid"] = mp_connection.get_mpid_from_formula(
                chem_formula=str(out_dict["formula"])
            )

        # Create a list out of a possible miller index string
        if "miller" in out_dict.keys() and isinstance(out_dict["miller"], str):
            out_dict["miller"] = [int(k) for k in list(out_dict["miller"])]

        # Check vdw and spin usage and correct to bool input:
        true_list = [
            "true",
            "True",
            "TRUE",
            ".TRUE.",
            ".true.",
            True,
            "yes",
            "Yes",
            "YES",
            ".YES.",
            ".yes.",
        ]
        vdw_keywords = [
            "dftd2",
            "dftd3",
            "dftd3-bj",
            "ts",
            "ts-hirshfeld",
            "mbd@rsc",
            "ddsc",
            "df",
            "optpbe",
            "optb88",
            "optb86b",
            "df2",
            "rvv10",
        ]
        if input_dict.get("use_vdw"):
            if input_dict.get("use_vdw") in true_list:
                out_dict["use_vdw"] = True
            elif input_dict.get("use_vdw") in vdw_keywords:
                out_dict["use_vdw"] = input_dict.get("use_vdw")
            else:
                out_dict["use_vdw"] = False
        if input_dict.get("use_spin"):
            if input_dict.get("use_spin") in true_list:
                out_dict["use_spin"] = True
            else:
                out_dict["use_spin"] = False

        return FWAction(mod_spec=[{"_set": {output_dict_name: out_dict}}])


# ============================================================================
# Initial handling of input dictionaries and parameters
# ============================================================================


def unbundle_input(
    inputs,
    keys=(
        "material_1",
        "material_2",
        "sg_params_1",
        "sg_params_2",
        "sg_filter_1",
        "sg_filter_2",
        "computational_params",
        "interface_params",
        "db_file",
        "high_level",
    ),
):
    """
    Read the input parameters and return the needed dictionaries for an
    interface made of two materials.

    Parameters
    ----------
    inputs : dict
        Dictionary containing the input parameters for a material. It is
        expected to have three keys, for the material, the computational, and
        the interfacial parameters, respectively.

    keys : tuple
        List containing the keys to be used to extract

    Returns
    -------
    material_1 : dict
        Parameters of the first material.

    material_2 : dict
        Parameters of the second material.

    comp_params : dict
        Computational parameters for the simulation.

    inter_params : dict
        Interface parameters to match the slabs.

    """

    # Get dicts containing material, computational and interfacial parameters
    # check if the set of keys in the inputs dictionary is a subset of the
    # required keys
    if not set(inputs.keys()).issubset(set(keys)):
        raise ValueError(
            "The inputs dictionary provided does not contain "
            "the correct keys:\n"
            "{}\n != \n{}".format(list(inputs.keys()), keys)
        )

    return inputs.values()


def material_from_mp(material_dict):
    """
    It reads the dictionary containing the input parameters for a material.
    It needs at least the key: `formula`, providing an MPID is also helpful.
    The corresponding structure from the MP database or the local structure list
    is returned along with the MPID.

    Parameters
    ----------
    material_dict : dict
        Dictionary containing the input parameters for running a tribological
        workflow, e.g. homogeneous and heterogeneous ones.

    Raises
    ------
    SystemExit
        If the parameters in the inputs dictionary are not set correctly.

    Returns
    -------
    struct : pymatgen.core.structure.Structure
        Contains the material structures.
    mpid : str or list of str
        Contains the MPID found for the corresponding structures
    functional : str
        The functional for the pseudopotential to be adopted.

    """

    # Collect the data from Material's Project API
    mp_connection = MPConnection()

    struct, mpid = mp_connection.get_low_energy_structure(
        chem_formula=material_dict.get("formula"),
        mp_id=material_dict.get("mpid"),
    )

    return struct, mpid


# ============================================================================
# Check the content of dictionary and make them coherent
# ============================================================================


def read_input_dict(input_dict, read_key):
    """
    Check the consistency of an input dictionary for a given read_key, compared
    with a dictionary containing a list of allowed keys and default values.

    Parameters
    ----------
    input_dict : dict
        Dictionary with input data to be analyzed.
    read_key : str
        Key to be read from defaults dictionary.

    Returns
    -------
    out_dict : dict
        Output dictionary after cleaning.

    """

    # Read the essential&additional list of keys and default values
    defaults = load_defaults()

    # Check if the read_key is known
    # if read_key not in list(defaults.keys()):
    #     raise ValueError(
    #         "There exist no default values for the selected key\n"
    #         "\t{}.\n Please choose from: {}".format(read_key, list(defaults.keys()))
    #     )
    # Extract the dictionary data
    out_dict = dict_consistency(input_dict, defaults[read_key])

    return out_dict


def dict_consistency(input_dict, defaults):
    """
    Check the consistency of an input dictionary. Dependency of read_input_dict.

    Parameters
    ----------
    input_dict : dict
        Input dictionary to be checked.
    defaults : dict
        Default values to be used when not specified differently from the user.

    Raises
    ------
    SystemExit
        If some essential keys are missing from input parameters or if some
        unknown values are used.

    Returns
    -------
    out_dict : dict
        Output dictionary after cleaning.

    """

    # file = open('output.txt', 'a')

    # Extract the input and parametric list of keys
    input_keys = list(input_dict.keys())
    essential_keys = []
    additional_keys = []
    for key, value in defaults.items():
        if value == "NO_DEFAULT":
            essential_keys.append(key)
        else:
            additional_keys.append(key)
    known_keys = essential_keys + additional_keys

    out_dict = {}

    # If essential keys are not a subset of inputs, something is missing
    if not set(essential_keys).issubset(set(input_keys)):
        raise SystemExit(
            "At least an essential input parameter is missing.\n"
            "\tInput parameters: {}\n"
            "\tEssential paremeters: {}\n".format(input_keys, essential_keys)
        )
        # file.write('At least an essential input parameter is missing. '
        #            'The following are required: ' + essential_keys)
        # raise SystemExit

    # If input keys are not a subset of known, at least one key is unknown
    if not set(input_keys).issubset(set(known_keys)):
        raise SystemExit(
            "At least an input parameter is not known.\n"
            "\tInput parameters: {} \n"
            "\tAllowed paremeters: {} \n".format(input_keys, known_keys)
        )
        # file.write('At least an input parameter is not known.\n'
        #            '\tInput parameters: ' + input_keys + '\n'
        #            '\tAllowed paremeters' + known_keys + '\n')
        # raise SystemExit

    # Read all input keys and use degaults for missing values
    for key in known_keys:
        out_dict[key] = input_dict.get(key, defaults.get(key, None))

    # file.close()

    return out_dict
