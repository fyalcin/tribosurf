#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:52:57 2020

@author: wolloch
"""
import numpy as np
import pickle
from atomate.utils.utils import env_chk
from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from monty.json import jsanitize
from operator import itemgetter
from pymatgen.core.interface import Interface
from pymatgen.core.structure import Structure

from htflow_utils.db_tools import VaspDB
from htflow_utils.vasp_tools import get_custom_vasp_relax_settings
from htflow_utils.workflows import dynamic_relax_swf
from triboflow.phys.high_symmetry import InterfaceSymmetryAnalyzer
from triboflow.phys.potential_energy_surface import (
    get_pes_generator_from_db,
)
from triboflow.utils.structure_manipulation import (
    clean_up_site_properties,
    get_interface_distance,
)


@explicit_serialize
class FT_ComputePES(FiretaskBase):
    """Compute the PES for a given interface, plot and save it.

    Uses the previously computed energies for the unique high-symmetry points
    and copies them to all the correct replica points. Replicates the points
    and fits the PES using radial basis functions. Output is saved in the
    database and if wanted also to files.

    :param interface_name: Name of the interface in the high-level database.
    :type interface_name: str

    :param functional: Which functional to use; has to be 'PBE' or 'SCAN'.
    :type functional: str

    :param external_pressure: External pressure on the interface in GPa.
    :type external_pressure: float

    :param file_output: Determines if results are written to disc.
    :type file_output: bool

    :param db_file: Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    :type db_file: str, optional

    :param high_level: Name of the high-level database.
    :type high_level: str, optional

    :return: None
    :rtype: NoneType
    """
    _fw_name = "Compute PES"
    required_params = [
        "interface_name",
        "functional",
        "file_output",
        "external_pressure",
    ]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        name = self.get("interface_name")
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)

        db_high = VaspDB(db_file=db_file, high_level=hl_db)

        data = db_high.find_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
        )
        if data:
            pes_data = data["PES"].get("all_energies", None)
            if pes_data:
                print("PES already computed for this interface. Skipping...")
                return FWAction()

        PG = get_pes_generator_from_db(
            interface_name=name,
            external_pressure=external_pressure,
            db_file=db_file,
            high_level=hl_db,
            functional=functional,
            pes_generator_kwargs={"fig_title": name},
        )

        # the data here might be too large to write to the DB, so if initial
        # write fails, we try bit by bit...
        # maybe should be replaced by gridFS instead?
        try:
            db_high.update_data(
                collection=functional + ".interface_data",
                fltr={"name": name, "external_pressure": external_pressure},
                new_values={
                    "$set": {
                        "PES.all_energies": jsanitize(PG.extended_energies),
                        "PES.pes_data": jsanitize(PG.PES_on_meshgrid),
                        "PES.image": PG.PES_as_bytes,
                        "PES.rbf": pickle.dumps(PG.rbf),
                        "corrugation": PG.corrugation,
                        "hsp@min": PG.hsp_min,
                        "hsp@max": PG.hsp_max,
                        "mep": jsanitize(PG.mep),
                        "shear_strength": jsanitize(PG.shear_strength),
                        "initial_strings": {
                            "x": PG.initial_string_x.tolist(),
                            "y": PG.initial_string_y.tolist(),
                            "d": PG.initial_string_d.tolist(),
                        },
                    }
                },
            )
        except:
            db_high.update_data(
                collection=functional + ".interface_data",
                fltr={"name": name, "external_pressure": external_pressure},
                new_values={
                    "$set": {
                        "corrugation": PG.corrugation,
                        "hsp@min": PG.hsp_min,
                        "hsp@max": PG.hsp_max,
                        "mep": jsanitize(PG.mep),
                        "shear_strength": jsanitize(PG.shear_strength),
                    }
                },
            )
            for k, v in {
                "all_energies": "extended_energies",
                "pes_data": "PES_on_meshgrid",
                "image": "PES_as_bytes",
                "rbf": "rbf",
            }.items():
                try:
                    if k == "rbf":
                        data = pickle.dumps(getattr(PG, v))
                    else:
                        data = jsanitize(getattr(PG, v))
                    db_high.update_data(
                        collection="PBE.interface_data",
                        fltr={
                            "name": name,
                            "external_pressure": external_pressure,
                        },
                        new_values={"$set": {"PES." + k: data}},
                    )
                except:
                    print(f"document {k} is too large to be put into DB.")
                    pass


@explicit_serialize
class FT_RetrievePESEnergies(FiretaskBase):
    """Retrieve the energies from the PES relaxations and update the db.

    Uses a tag together with the labels of the high-symmetry points saved
    in the high level database to retrieve the correct energies for each
    lateral shift of the interface. Sort the shifts by energies and save both
    the configuration with the lowest and the highest energy in the high level
    database. Also save the list of shifts and energies with corresponding
    labels there.

    :param interface_name: Name of the interface in the high-level database.
    :type interface_name: str

    :param functional: Which functional to use; has to be 'PBE' or 'SCAN'.
    :type functional: str

    :param tag: Unique tag to identify the calculations.
    :type tag: str

    :param external_pressure: External pressure on the interface in GPa.
    :type external_pressure: float

    :param db_file: Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    :type db_file: str, optional

    :return: None
    :rtype: NoneType
    """
    _fw_name = "Retrieve PES energies"
    required_params = [
        "interface_name",
        "functional",
        "tag",
        "external_pressure",
    ]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        name = self.get("interface_name")
        functional = self.get("functional")
        tag = self.get("tag")
        external_pressure = self.get("external_pressure")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)

        db_high = VaspDB(db_file=db_file, high_level=hl_db)

        interface_dict = db_high.find_data(
            collection=f"{functional}.interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
        )
        lateral_shifts = interface_dict["PES"]["high_symmetry_points"][
            "unique_shifts"
        ]
        group_assignments = interface_dict["PES"]["high_symmetry_points"][
            "group_assignments"
        ]

        db = VaspDB(db_file=db_file)
        ref_struct = Interface.from_dict(
            interface_dict["unrelaxed_structure"]
        )
        area = np.linalg.norm(
            np.cross(
                ref_struct.lattice.matrix[0], ref_struct.lattice.matrix[1]
            )
        )

        energy_list = []
        energy_dict = {}
        calc_output = {}
        for s in lateral_shifts.keys():
            label = tag + "_" + s
            vasp_calc = db.find_data(
                collection="tasks", fltr={"task_label": label}
            )
            struct = vasp_calc["output"]["structure"]
            energy = vasp_calc["output"]["energy"]
            energy *= 16.02176565 / area
            energy_list.append([s, group_assignments[s], energy])
            energy_dict[s] = energy
            calc_output[s] = {
                "energy": energy,
                "relaxed_struct": struct,
                "task_id": vasp_calc["_id"],
            }

        sorted_energy_list = sorted(energy_list, key=itemgetter(2))

        min_stacking = sorted_energy_list[0][0]
        max_stacking = sorted_energy_list[-1][0]
        calc_min = db.find_data(
            collection="tasks",
            fltr={"task_label": tag + "_" + min_stacking},
        )
        calc_max = db.find_data(
            collection="tasks",
            fltr={"task_label": tag + "_" + max_stacking},
        )
        struct_min_dict = calc_min["output"]["structure"]
        struct_max_dict = calc_max["output"]["structure"]

        struct_min = Structure.from_dict(struct_min_dict)
        struct_max = Structure.from_dict(struct_max_dict)

        for index in range(ref_struct.num_sites):
            struct_min[index].properties = ref_struct[index].properties
            struct_max[index].properties = ref_struct[index].properties

        struct_min_dict = struct_min.as_dict()
        struct_max_dict = struct_max.as_dict()
        struct_min_dict.update(
            {"interface_properties": ref_struct.interface_properties}
        )
        struct_max_dict.update(
            {"interface_properties": ref_struct.interface_properties}
        )

        interface_min = Interface.from_dict(struct_min_dict)
        interface_max = Interface.from_dict(struct_max_dict)

        # add site properties from ref_struct to get back interface_labels:
        for k, v in ref_struct.site_properties.items():
            struct_min.add_site_property(k, v)
            struct_max.add_site_property(k, v)

        inter_dist_min = get_interface_distance(struct_min)
        inter_dist_max = get_interface_distance(struct_max)

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        db_high.update_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
            new_values={
                "$set": {
                    "relaxed_structure@min": interface_min.as_dict(),
                    "relaxed_structure@max": interface_max.as_dict(),
                    "interface_distance@min": inter_dist_min,
                    "interface_distance@max": inter_dist_max,
                    "PES.calculations": calc_output,
                    "PES.high_symmetry_points.energies_dict": energy_dict,
                }
            },
        )


@explicit_serialize
class FT_FindHighSymmPoints(FiretaskBase):
    """Compute high symmetry points for the top and bottom slab and the interface.

    Finds the high symmetry points of the top side of the bottom slab and the
    bottom side of the top slab. This is done twice, once omitting duplicates,
    and once allowing them. It is made sure that the results are cartesian
    coordinates that lie inside the unit cell. The lists are combined so that
    every combination of unique points for the interface is present. The PES
    section of the high level database is updated with the results and the
    fw_spec is updated with the lateral shifts needed for the PES relaxations
    as well.

    :param interface: Interface object for which the PES is to be constructed
    :type interface: pymatgen.core.interface.Interface

    :param functional: Which functional to use; has to be 'PBE' or 'SCAN'.
    :type functional: str

    :param interface_name: Name of the interface in the high-level database.
    :type interface_name: str

    :param external_pressure: External pressure in GPa.
    :type external_pressure: float

    :param db_file: Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    :type db_file: str, optional

    :return: FWActions that updates the fw_spec with lateral shifts.
    :rtype: FWAction
    """

    required_params = [
        "interface",
        "functional",
        "interface_name",
        "external_pressure",
    ]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        interface = self.get("interface")
        name = self.get("interface_name")
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        interface_data = db_high.find_data(
            f"{functional}.interface_data",
            {"name": name, "external_pressure": external_pressure},
        )
        interface_params = interface_data["interface_params"]

        ISA = InterfaceSymmetryAnalyzer(interface)
        hsp_dict = ISA.get_high_symmetry_info()
        interfaces = ISA.get_all_high_symmetry_interfaces(
            interface_distance=interface_params.get("interface_distance")
        )

        db_high.update_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
            new_values={
                "$set": {
                    "PES.high_symmetry_points": {
                        "bottom_unique": hsp_dict[
                            "bottom_high_symm_points_unique"
                        ],
                        "bottom_all": hsp_dict[
                            "bottom_high_symm_points_all"
                        ],
                        "top_unique": hsp_dict[
                            "top_high_symm_points_unique"
                        ],
                        "top_all": hsp_dict["top_high_symm_points_all"],
                        "unique_shifts": hsp_dict["unique_shifts"],
                        "all_shifts": hsp_dict["all_shifts"],
                        "group_assignments": hsp_dict["group_assignments"],
                    }
                }
            },
            upsert=True,
        )
        fw_spec["high_symm_interfaces"] = interfaces

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartPESCalcs(FiretaskBase):
    """Start z-relaxations for different lateral positions of an interface.

    Take a list of lateral shifts from the fw_spec and start relaxations
    for each one of them as parallel detours.


    :param interface_name: Name of the interface in the high-level database.
    :type interface_name: str

    :param comp_parameters: Computational parameters to be passed to the vasp input
        file generation.
    :type comp_parameters: dict

    :param tag: Unique tag to identify the calculations.
    :type tag: str

    :param external_pressure: External pressure on the interface in GPa.
    :type external_pressure: float

    :param db_file: Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    :type db_file: str, optional

    :param prerelax: Whether to perform a prerelaxation using a network potential
        before starting a DFT relaxation. Defaults to True.
    :type prerelax: bool, optional

    :param prerelax_calculator: Which network potential to use for the prerelaxation.
        Defaults to 'm3gnet'.
    :type prerelax_calculator: str, optional

    :param prerelax_kwargs: Keyword arguments to be passed to the ASE calculator for
        the prerelaxation.
    :type prerelax_kwargs: dict, optional

    :return: FWActions that produce a detour workflow with relaxations for the PES.
    :rtype: FWAction
    """
    _fw_name = "Start PES calculations"
    required_params = [
        "interface_name",
        "comp_parameters",
        "tag",
        "external_pressure",
    ]
    optional_params = [
        "prerelax",
        "prerelax_calculator",
        "prerelax_kwargs",
        "db_file",
    ]

    def run_task(self, fw_spec):
        interface_name = self.get("interface_name")
        comp_params = self.get("comp_parameters")
        tag = self.get("tag")
        external_pressure = self.get("external_pressure")

        prerelax = self.get("prerelax", True)
        prerelax_calculator = self.get("prerelax_calculator", "m3gnet")
        prerelax_kwargs = self.get("prerelax_kwargs", {})
        db_file = self.get("db_file", "auto")

        interfaces = fw_spec.get("high_symm_interfaces")
        if not interfaces:
            raise SystemExit(
                "High symmetry interface list not found in the fw_spec./n"
                "Please check your Firework for errors!"
            )

        inputs = []
        for interface in interfaces:
            group_name = interface.interface_properties[
                "high_symmetry_info"
            ]["group_name"]
            label = tag + "_" + group_name
            clean_struct = clean_up_site_properties(interface)

            vis = get_custom_vasp_relax_settings(
                structure=clean_struct,
                comp_parameters=comp_params,
                relax_type="interface_z_relax",
                apply_pressure=external_pressure,
            )
            inputs.append([clean_struct, vis, label])

        wf_name = "PES relaxations for: " + interface_name
        WF = dynamic_relax_swf(
            inputs_list=inputs,
            wf_name=wf_name,
            add_static=True,
            prerelax_system=prerelax,
            prerelax_calculator=prerelax_calculator,
            prerelax_kwargs=prerelax_kwargs,
            db_file=db_file,
        )

        return FWAction(detours=WF)
