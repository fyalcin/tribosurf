#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:16:55 2021

Collection of Firetasks to start workflows concerning slabs structure.
They can be used to run existing workflows as subworkflows within other ones.

The module contains the following Firetasks:

    - FT_SlabOptThick 
    Starts a subworkflow to perform a convergence process to find the optimal 
    thickness of a slab structure. First step to be called within a workflow.

    - FT_StartThickConvo
    Start a subworkflow to select a desired convergence criterion to calculate
    the slab thickness. Implemented criteria: surface energy.

    - FT_EndThickConvo
    Check the results and call recursively the slab thickness workflow it the
    convergence has not been achieved yet.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Gabriele Losi"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 22nd, 2021"


import os
from monty.json import jsanitize
from hitmen_utils.db_tools import VaspDB
import numpy as np
from pymatgen.core.structure import Structure
from fireworks import explicit_serialize, FiretaskBase, FWAction

from triboflow.utils.database import (
    Navigator,
    StructureNavigator,
    get_low_and_high_db_names,
)
from triboflow.utils.utils import (
    read_runtask_params,
    read_default_params,
    get_one_info_from_dict,
    retrieve_from_db,
)
from triboflow.utils.structure_manipulation import slab_from_structure
from triboflow.utils.errors import SlabOptThickError
from triboflow.phys.shaper import Shaper
from hitmen_utils.misc_tools import check_input

currentdir = os.path.dirname(__file__)


@explicit_serialize
class FT_SlabOptThick(FiretaskBase):
    """
    Start a subworkflow as a detour to calculate the optimal thickness for a
    slab generated from a supplied bulk structure and a certain orientation
    (miller index). The thickness can be converted either by evaluating how the
    surface energy or the lattice parameters (not implemented) change with the
    number of atomic planes.

    Parameters
    ----------
    mp_id : str
        MP-id of the structure from the MP database.

    miller : list of int, or str
        Miller indexes (h, k, l) to select the slab orientation.

    functional : str
        Functional for the pseudopotential to be adopted.

    db_file : str or None
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.

    low_level : str or None, optional
        Name of the table of the "low level" database, saved in db_file. The
        intermediate calculations and raw data during will be saved here. If
        nothing is passed the Firework database is used. The default is None.

    high_level : str, or None, optional
        Name of the table of the "high level" database, saved in db_file.
        The slab optimal thickness will be saved here. The slab energy and
        surface energy will be saved too if conv_kind is 'surfene'.
        The default is None.

    conv_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene',
        'alat'. The latter is not implemented yet. The default is 'surfene'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to `get_custom_vasp_relax_settings`. The default is 'slab_pos_relax'.

    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point. In
        case of low parallelization this gives the value of atomic layers for
        the only slab which is created and relaxed at each step after the first.
        The default is 4.

    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If convergence is
        not reached this value is considered the optimal one. The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step during
        the iterative procedure. The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.

    ext_index : int, optional
        Use the ext_index element from SlabGenerator.get_slabs as a slab.
        The default is 0.

    conv_thr : float, optional
        Threshold for the convergence process. The default is 0.025.

    bulk_entry : str or list or None, optional
        Name of the custom bulk dictionary, to be retrieved from the high level
        database and to be used to build the slabs. Bulks are identified by
        mp_id and functional but there might be different structures of the
        same material. The default is "structure_fromMP".

    cluster_params : dict, optional
        Optional params to print data and/or interact with clusters. The default is {}.

    override : bool, optional
        Decide if the dft simulation should be done in any case, despite the
        possible presence of previous results.

    add_static : bool, optional
       Selects if a static calculation is done after the relaxation. This
       is useful for accurate total energies.

    """

    _fw_name = "Start a subworkflow to converge slab thickness"

    required_params = ["mp_id", "miller", "functional"]
    optional_params = [
        "db_file",
        "low_level",
        "high_level",
        "conv_kind",
        "relax_type",
        "thick_min",
        "thick_max",
        "thick_incr",
        "vacuum",
        "in_unit_planes",
        "ext_index",
        "conv_thr",
        "parallelization",
        "bulk_entry",
        "cluster_params",
        "override",
        "add_static",
    ]

    def run_task(self, fw_spec):
        """Run the Firetask."""
        if type(self.get("miller")) == str:
            miller = [int(k) for k in list(self.get("miller"))]
        else:
            miller = self.get("miller")
        nav = StructureNavigator(db_file="auto", high_level=True)
        slab_data = nav.get_slab_from_db(
            self.get("mp_id"), self.get("functional"), miller
        )
        for key in self.optional_params:
            if slab_data["slab_parameters"].get(key):
                self[key] = slab_data["slab_parameters"].get(key)

        if not self.get("conv_thr"):
            self["conv_thr"] = slab_data.get("comp_parameters").get(
                "surfene_thr"
            )
            print(self["conv_thr"])

        # Define the json file containing default values and read parameters
        dfl = currentdir + "/../defaults.json"
        p = read_runtask_params(
            self,
            fw_spec,
            self.required_params,
            self.optional_params,
            default_file=dfl,
            default_key="SlabOptThick",
        )

        # Check if a convergence calculation of the slab thickness is present
        is_done, comp_params = self.is_data(p, dfl)

        import pprint

        pprint.pprint(p)

        # Start a subworkflow to converge the thickness if not already done
        if is_done and not p["override"]:
            return FWAction(update_spec=fw_spec)  # Continue the Workflow

        else:
            # Create the subworkflow to run a convergence detour
            bulk = self.get_bulk(p)
            wf = self.select_slabthick_conv(
                structure=bulk,
                fw_spec=fw_spec,
                comp_params=comp_params,
                dfl=dfl,
                p=p,
            )

            return FWAction(detours=wf, update_spec=fw_spec)

    def get_bulk(self, p):
        """
        Get the bulk structure from the database, to generate the slabs.

        Parameters
        ----------
        p : dict
            Input parameters of the Firetask.

        Returns
        -------
        bulk : pymatgen.core.structure.Structure
            pymatgen bulk structure
        """

        field, bulk = retrieve_from_db(
            p["mp_id"],
            collection=p["functional"] + ".bulk_data",
            db_file=p["db_file"],
            high_level_db=True,
            entry=p["bulk_entry"],
            is_slab=False,
            pymatgen_obj=True,
        )
        return bulk

    def is_data(self, p, dfl):
        """
        Query the database, download the slab dictionary and check if an entry
        for the optimal thickness already exists.

        Parameters
        ----------
        p : dict
            Input parameters of the Firetask.

        Returns
        -------
        is_data : bool
            True if a key named 'opt_thickness' is found in `entry`

        comp_params : dict
            Computational parameters to simulate the slab.

        """

        # Retrieve the slab from the database
        field, _ = retrieve_from_db(
            mp_id=p["mp_id"],
            db_file=p["db_file"],
            high_level_db=True,
            miller=p["miller"],
            collection=p["functional"] + ".slab_data",
            pymatgen_obj=False,
        )

        if field is not None:
            # Check if an optimal thickness has been already calculated
            is_done = field.get("opt_thickness", None)

            # Define the computational parameters of the slab
            comp_params = field.get("comp_parameters", {})

        else:
            is_done = False
            comp_params = {}

        return bool(is_done), comp_params

    def select_slabthick_conv(self, structure, fw_spec, comp_params, dfl, p):
        """
        Select the desired subworkflow from the SlabWFs class, to converge the
        slab thickness either by evaluating the surface energy or the lattice
        parameter.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            Bulk structure to construct the slabs from.

        fw_spec : dict
            Spec of the dictionary, contains information to bootstrap job.
            It is not so relevant in the present workflow, however it might be
            important in future workflows containing the SlabThickOpt Firetask
            to keep information through this convergence process.

        comp_params : dict
            Computational parameters of the slab.

        dfl : str
            Path to the JSON file containing default values for the workflows.

        p : dict
            Input parameters of the Firetask.

        Returns
        -------
        wf : Firework.Workflow
            Workflow object to start a convergence process as a detour.

        """

        from triboflow.workflows.slabs_wfs import SlabWF

        # Check for default values of comp_params and cluster_params
        comp_params = read_default_params(dfl, "comp_params", comp_params)
        cluster_params = read_default_params(
            dfl, "cluster_params", p["cluster_params"]
        )

        # Select the convergence function based on conv_kind
        if p["conv_kind"] == "surfene":
            generate_wf = SlabWF.conv_slabthick_surfene
        elif p["conv_kind"] == "alat":
            generate_wf = SlabWF.conv_slabthick_alat
        else:
            raise SlabOptThickError(
                "Wrong input argument for conv_kind. "
                "Allowed options: 'surfene', 'alat'"
            )

        if structure is None:
            raise SlabOptThickError(
                "Bulk has not been found in database:\n"
                "db_file : {}\ndatabase : {}\ncollection : {}\n"
                "entry : {}".format(
                    p["db_file"],
                    p["high_level"],
                    p["functional"] + ".bulk_data",
                    p["bulk_entry"],
                )
            )

        # Generate the workflow
        wf = generate_wf(
            structure=structure,
            mp_id=p["mp_id"],
            spec=fw_spec,
            miller=p["miller"],
            functional=p["functional"],
            comp_params=comp_params,
            db_file=p["db_file"],
            low_level=p["low_level"],
            high_level=p["high_level"],
            relax_type=p["relax_type"],
            thick_min=p["thick_min"],
            thick_max=p["thick_max"],
            thick_incr=p["thick_incr"],
            vacuum=p["vacuum"],
            in_unit_planes=p["in_unit_planes"],
            ext_index=p["ext_index"],
            conv_thr=p["conv_thr"],
            parallelization=p["parallelization"],
            cluster_params=cluster_params,
            override=p["override"],
            add_static=p["add_static"],
        )

        return wf


@explicit_serialize
class FT_StartThickConvo(FiretaskBase):
    """
    It starts a subworkflow as a detour to converge the thickness of a slab.
    The thickness can be converged either by evaluating how the surface energy
    or the lattice parameters changes with the number of atomic planes.
    It is the first element of the SlabWF.conv_slabthick_surfene workflow.

    (At the moment only the surface energy convergence is implemented)

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Bulk structure to build the slabs.

    mp_id : str
        MP-id of the structure from the MP database.

    miller : list of int, or str
        Miller indexes (h, k, l) to select the slab orientation.

    functional : str
        Functional for the pseudopotential to be adopted.

    db_file : str or None
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.

    low_level : str or None, optional
        Name of the table of the "low level" database, saved in db_file. The
        intermediate calculations and raw data during will be saved here. If
        nothing is passed the Firework database is used. The default is None.

    high_level : str, optional
        Name of the table of the "high level" database, saved in db_file.
        The slab optimal thickness will be saved here. The slab energy and
        surface energy will be saved too if conv_kind is 'surfene'.
        The default is 'triboflow'.

    conv_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene',
        'alat'. The latter is not implemented yet. The default is 'surfene'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.

    comp_params : dict, optional
        Computational parameters for the VASP simulations. If not set, default
        parameters will be used instead. The default is {}.

    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point. In
        case of low parallelization this gives the value of atomic layers for
        the only slab which is created and relaxed at each step after the first.
        The default is 4.

    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If convergence is
        not reached this value is considered the optimal one. The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step during
        the iterative procedure. The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.

    ext_index : int, optional
        Use the ext_index element from SlabGenerator.get_slabs as a slab.
        The default is 0.

    cluster_params : dict, optional
        Dictionary containing cluster-related options to run efficiently the
        VASP simulations on a cluster. The default is {}.

    """

    _fw_name = "Start the slab thickness convergence"
    required_params = ["structure", "mp_id", "miller", "functional"]
    optional_params = [
        "db_file",
        "low_level",
        "high_level",
        "conv_kind",
        "relax_type",
        "comp_params",
        "thick_min",
        "thick_max",
        "thick_incr",
        "vacuum",
        "in_unit_planes",
        "ext_index",
        "parallelization",
        "recursion",
        "cluster_params",
        "override",
        "add_static",
    ]

    def run_task(self, fw_spec):
        """Run the Firetask."""

        # Define the json file containing default values and read parameters
        dfl = currentdir + "/../defaults.json"
        p = read_runtask_params(
            self,
            fw_spec,
            self.required_params,
            self.optional_params,
            default_file=dfl,
            default_key="StartThickConvo",
        )

        # Select the convergence of interest
        wf = self.select_conv(p)

        return FWAction(detours=wf, update_spec=fw_spec)

    def select_conv(self, p):
        from triboflow.workflows.surfene_wfs import SurfEneWF

        if p["conv_kind"] == "surfene":
            low_level_name, high_level_name = get_low_and_high_db_names(p)
            wf = SurfEneWF.conv_surface_energy(
                structure=p["structure"],
                mp_id=p["mp_id"],
                miller=p["miller"],
                functional=p["functional"],
                db_file=p["db_file"],
                low_level=low_level_name,
                high_level=high_level_name,
                relax_type=p["relax_type"],
                comp_params=p["comp_params"],
                thick_min=p["thick_min"],
                thick_max=p["thick_max"],
                thick_incr=p["thick_incr"],
                vacuum=p["vacuum"],
                in_unit_planes=p["in_unit_planes"],
                ext_index=p["ext_index"],
                parallelization=p["parallelization"],
                recursion=p["recursion"],
                cluster_params=p["cluster_params"],
                override=p["override"],
                add_static=p["add_static"],
            )
            return wf

        else:
            raise SystemExit(
                "Lattice parameter convergence not yet implemented"
            )


@explicit_serialize
class FT_EndThickConvo(FiretaskBase):
    """
    Firetask description...

    """

    required_params = ["structure", "mp_id", "miller"]
    optional_params = [
        "db_file",
        "low_level",
        "high_level",
        "spec",
        "functional",
        "conv_kind",
        "relax_type",
        "comp_params",
        "thick_min",
        "thick_max",
        "thick_incr",
        "vacuum",
        "in_unit_planes",
        "ext_index",
        "conv_thr",
        "parallelization",
        "cluster_params",
        "recursion",
        "override",
        "add_static",
    ]

    def run_task(self, fw_spec):
        """Run the Firetask."""

        # Define the json file containing default values and read parameters
        dfl = currentdir + "/../defaults.json"
        p = read_runtask_params(
            self,
            fw_spec,
            self.required_params,
            self.optional_params,
            default_file=dfl,
            default_key="EndThickConvo",
        )
        case = self._get_case(p)

        # Retrieve the surface energy
        data, index = self.get_data(case, p)

        # Analyze the data, if convergence is reached data is stored in DB
        stop_convergence = self.analyze_data(data, index, case, p)

        # Rerun recursively the surface energy subworkflow within detours
        if not stop_convergence:
            wf = self.call_recursion(fw_spec, p)

            return FWAction(detours=wf, update_spec=fw_spec)

        else:
            return FWAction(update_spec=fw_spec)

    def _get_case(self, p):
        """
        Define the dictionary key to be read from 'conv_kind'.

        """

        # Select the keys to be used when reading the dictionary
        if p["conv_kind"] == "surfene":
            case = "surface_energy"
        elif p["conv_kind"] == "alat":
            case = "lattice"
        else:
            raise SlabOptThickError(
                "Wrong argument: 'conv_kind'. Allowed "
                "values: 'surfene', 'alat'"
            )
        return case

    def get_data(self, case, p):
        """
        Extract the surface energies from the high level DB.

        case : str
            Dictionary key to identify the type of data to be read from
            'calc_output' in the nested entry of the pymongo field.

        p : dict
            Dictionary with the input parameters of the Firetask.

        Returns
        -------
        data : list of floats
            Contains all the `case` data calculated for various thicknesses.

        index : list of index
            Contains the index referring to the different thicknesses.

        """

        # Call the navigator for retrieving the thickness dict out of DB
        nav = Navigator(db_file=p["db_file"], high_level=False)
        thickness_dict = nav.find_data(
            p["functional"] + ".slab_data",
            {"mpid": p["mp_id"], "miller": p["miller"]},
        )["thickness"]

        # Get the indexes with the thickness and the desired data
        index = []
        data = []

        # Define the thicknesses that should be taken into account
        thks = []
        if p["parallelization"] == "high":
            for key in thickness_dict.keys():
                if key.startswith("data_"):
                    i = int(key.split("_")[-1])
                    if i != 0:
                        thks.append(i)
            if not p["thick_min"] in thks or not p["thick_max"] in thks:
                raise SlabOptThickError(
                    "Min or Max values not present among "
                    "the results. Something weird happening"
                )
        else:
            thks = [p["thick_min"], p["thick_max"]]

        # Read allowed data for this turn
        for key, item in thickness_dict.items():
            if key.startswith("data_"):
                # Save the data of interest
                i = int(key.split("_")[-1])
                if i in thks:
                    if not "output" in item.keys():
                        raise SlabOptThickError(
                            "The calculation with {} layers "
                            "has no output, something "
                            "went wrong".format(i)
                        )
                    index.append(i)
                    print(item["output"].keys())
                    data.append(item["output"][case])

        sorted_index = np.array(index).argsort()
        index = np.array(index)[sorted_index]
        data = np.array(data)[sorted_index]

        # Check for consistency with the value from the max thickness
        dmax = thickness_dict["data_" + str(p["thick_max"])]["output"][case]
        if dmax != data[-1]:
            raise SlabOptThickError("An unexpected error occurred")

        return data, index

    def analyze_data(self, data, index, case, p):
        """
        Analyze the data to understand if the convergence has been achieved.
        Put it in the `high_level` db at the end of the process.

        Parameters
        ----------
        data : list of floats
            List containing the values to be checked for convergence.

        index : list of int
            List containing the thicknesses of the different slabs.

        case : str
            Dictionary key to identify the type of data and store them.

        p : dict
            Dictionary with the input parameters of the Firetask.

        Returns
        -------
        stop_convergence : bool
            Decide the fate of the calculation.

        """

        # Calculate the relative error to the last element
        error_to_last = np.abs((data - data[-1]) / data[-1])

        # Evaluate what is the lower converged data
        i = np.argwhere(error_to_last <= abs(p["conv_thr"]))
        index_converged = index[i]
        index_converged = index_converged.flatten()

        # If a low parallelization is selected, three calculations are done
        # at the beginning (bulk, min, max), and if convergence is not
        # achieved, then a calculation is done one by one until convergence
        # is not attained or you go above max
        if p["parallelization"] in [None, "low"]:
            # If length>1, then a material has converged other than max
            if len(index_converged) > 1:
                self.store_to_db(index_converged[0], p)
                stop_convergence = True

            # If length=1, material is not converged, but if the next step
            # brings you above the max thickness, then convergence is
            # assumed to be reached at thick_max
            elif p["thick_min"] + p["thick_incr"] >= p["thick_max"]:
                self.store_to_db(index_converged[0], p)
                stop_convergence = True

            # If length=1 and you are far from thick max, recursion is done
            else:
                stop_convergence = False

        # If high level parallelization
        elif p["parallelization"] == "high":
            self.store_to_db(index_converged[0], p)
            stop_convergence = True

        return stop_convergence

    def store_to_db(self, index, p):
        """
        Store the surface energy or the latice parameter that have been
        obtained in a more accessible point within the database field.
        If the same convergence calculation is called again, it will stop
        when it finds that the parameters is already converged.

        Parameters
        ----------
        index : int
            Index of the thickness indicating the optimal slab.

        p : dict
            Dictionary containing the FT arguments.

        """

        # Start the navigator to extract the data from the low level database
        nav_low = Navigator(db_file=p["db_file"], high_level=False)
        low_dict = nav_low.find_data(
            collection=p["functional"] + ".slab_data",
            fltr={"mpid": p["mp_id"], "miller": p["miller"]},
        )

        # Start the navigator to store the data to the high level database
        nav_high = Navigator(db_file=p["db_file"], high_level=True)
        high_dict = nav_high.find_data(
            collection=p["functional"] + ".slab_data",
            fltr={"mpid": p["mp_id"], "miller": p["miller"]},
        )

        # Extract the data to be saved in the database
        thickness_dict = get_one_info_from_dict(low_dict, ["thickness"])
        out_struct_dict = thickness_dict["data_" + str(index)]["output"][
            "structure"
        ]

        # Convert out_struct_dict to slab (Issue in Atomate OptimizeFW, every
        # structure that is simulated in that way is saved in tasks collection
        # as a Structure)
        output_slab = slab_from_structure(
            p["miller"], Structure.from_dict(out_struct_dict)
        )

        opt_thickness = {
            "layers": int(index),
            "angstroms": Shaper.get_proj_height(output_slab, "slab"),
        }

        # Create an array containing the thickness vs surface energy info
        thick_array = []
        surfene_array = []
        for k in thickness_dict.keys():
            thick = int(k.split("_")[-1])
            if thick != 0:
                thick_array.append(thick)
                surfene_array.append(
                    thickness_dict[k]["output"]["surface_energy"]
                )

        array = np.column_stack((thick_array, surfene_array))
        thickness_dict["thick_surfene_array"] = array[np.argsort(array[:, 0])]
        opt_surfen = thickness_dict[f"data_{index}"]["output"][
            "surface_energy"
        ]
        opt_surfen_dict = {
            "eV/A^2": opt_surfen * 6.241509e-2,
            "J/m^2": opt_surfen,
        }

        # Prepare the dictionary for the update
        if high_dict is None:
            store = {
                "formula": output_slab.composition.reduced_formula,
                "mpid": p["mp_id"],
                "miller": p["miller"],
                "thickness": thickness_dict,
                "opt_thickness": opt_thickness,
                "relaxed_slab": output_slab.as_dict(),
                "surface_energy": opt_surfen_dict,
            }
        else:
            store = {
                "thickness": thickness_dict,
                "opt_thickness": opt_thickness,
                "relaxed_slab": output_slab.as_dict(),
                "surface_energy": opt_surfen_dict,
            }
        store = jsanitize(store)

        # Update data
        nav_high.update_data(
            collection=p["functional"] + ".slab_data",
            fltr={"mpid": p["mp_id"], "miller": p["miller"]},
            new_values={"$set": store},
            upsert=True,
        )

    def call_recursion(self, fw_spec, p):
        """
        Call the convergence workflow on the slab optimal thickness in a
        recursive way.

        """

        from triboflow.workflows.slabs_wfs import SlabWF

        # Select the correct function to call the workflow
        if p["conv_kind"] == "surfene":
            generate_wf = SlabWF.conv_slabthick_surfene
        else:
            pass

        # Generate the workflow for the detour
        wf = generate_wf(
            structure=p["structure"],
            mp_id=p["mp_id"],
            miller=p["miller"],
            functional=p["functional"],
            comp_params=p["comp_params"],
            spec=fw_spec,
            db_file=p["db_file"],
            low_level=p["low_level"],
            high_level=p["high_level"],
            relax_type=p["relax_type"],
            thick_min=p["thick_min"] + p["thick_incr"],
            thick_max=p["thick_max"],
            thick_incr=p["thick_incr"],
            vacuum=p["vacuum"],
            in_unit_planes=p["in_unit_planes"],
            ext_index=p["ext_index"],
            conv_thr=p["conv_thr"],
            parallelization=p["parallelization"],
            recursion=p["recursion"] + 1,
            override=p["override"],
            cluster_params=p["cluster_params"],
            add_static=p["add_static"],
        )

        return wf


@explicit_serialize
class GetSurfaceEnergiesFromUids(FiretaskBase):
    _fw_name = "Query the surface energies with the list of uids provided"
    required_params = ["uids"]
    optional_params = [
        "db_file",
        "high_level",
        "fake_calc",
        "surfen_coll"
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        nav = VaspDB(db_file=inp["db_file"], high_level=inp["high_level"])
        results = list(nav.find_many_data(inp["surfen_coll"], {"uid": {"$in": inp["uids"]}}))
        uid_surfen_dict = {r["uid"]: r["surface_energy"] for r in results}
        return FWAction(update_spec={"uid_surfen_dict": uid_surfen_dict})