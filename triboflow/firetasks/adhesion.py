#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 15:36:52 2021

@author: wolloch
"""
import numpy as np
from atomate.utils.utils import env_chk
from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab

from hitmen_utils.db_tools import VaspDB
from hitmen_utils.misc_tools import make_calculation_hash
from hitmen_utils.vasp_tools import get_custom_vasp_relax_settings
from hitmen_utils.workflows import dynamic_relax_swf
from triboflow.utils.structure_manipulation import (
    slab_from_structure,
)


@explicit_serialize
class FT_RetrieveMatchedSlabs(FiretaskBase):
    """Retrieve relaxed matched slabs and save them in the high_level database.

    Get the relaxed aligned top and bottom slabs from the low level database
    and save them in the interface_data collection of the high level database.

    :param mp_id_1: MaterialsProject ID number for the first material
    :type mp_id_1: str

    :param mp_id_2: MaterialsProject ID number for the second material
    :type mp_id_2: str

    :param functional: Functional with which the workflow is run. PBE or SCAN.
    :type functional: str

    :param external_pressure: External pressure in GPa.
    :type external_pressure: float

    :param db_file: Full path of the db.json file to be used. The default is to use env_chk to find the file.
    :type db_file: str, optional

    :param top_out_name: Name the relaxed top slab will have in the high-level database. Defaults to 'top_aligned_relaxed'
    :type top_out_name: str, optional

    :param bottom_out_name: Name the relaxed bottom slab will have in the high-level database. Defaults to 'bottom_aligned_relaxed'
    :type bottom_out_name: str, optional

    :param high_level: Name of the high_level database to use. Defaults to 'True', in which case it is read from the db.json file.
    :type high_level: str or True, optional

    :return: FWAction that updates the spec.
    :rtype: FWAction
    """
    _fw_name = "Retrieve relaxed matched slabs"
    required_params = [
        "mp_id_1",
        "mp_id_2",
        "functional",
        "external_pressure",
    ]
    optional_params = [
        "db_file",
        "top_out_name",
        "bottom_out_name",
        "high_level",
    ]

    def run_task(self, fw_spec):
        mp_id_1 = self.get("mp_id_1")
        mp_id_2 = self.get("mp_id_2")
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        top_out_name = self.get("top_out_name", "top_aligned_relaxed")
        bot_out_name = self.get("bottom_out_name", "bottom_aligned_relaxed")
        hl_db = self.get("high_level", True)

        name = fw_spec.get("interface_name")
        input_list = fw_spec.get("relaxation_inputs")

        if input_list:
            db = VaspDB(db_file)
            for i in input_list:
                label = i[-1]
                miller = i[0].miller_index
                calc = db.find_data(
                    collection="tasks", fltr={"task_label": label}
                )
                out_struct = calc["output"]["structure"]
                slab = slab_from_structure(
                    miller, Structure.from_dict(out_struct)
                )
                if label.startswith("top"):
                    out_name = top_out_name
                else:
                    out_name = bot_out_name
                db_high = VaspDB(db_file, high_level=hl_db)
                db_high.update_data(
                    collection=functional + ".interface_data",
                    fltr={
                        "name": name,
                        "external_pressure": external_pressure,
                    },
                    new_values={"$set": {out_name: slab.as_dict()}},
                )

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_RelaxMatchedSlabs(FiretaskBase):
    """Start the relaxation of the matched slabs.

    Get the aligned top and bottom slabs from the high level database and
    start their relaxation runs.

    :param mp_id_1: MaterialsProject ID number for the first material
    :type mp_id_1: str

    :param mp_id_2: MaterialsProject ID number for the second material
    :type mp_id_2: str

    :param functional: Functional with which the workflow is run. PBE or SCAN.
    :type functional: str

    :param external_pressure: External pressure in GPa.
    :type external_pressure: float

    :param db_file: Full path of the db.json file to be used. The default is to use env_chk to find the file.
    :type db_file: str, optional

    :param top_in_name: Name the unrelaxed top slab has in the high-level database. Defaults to 'top_aligned'
    :type top_in_name: str, optional

    :param bottom_in_name: Name the unrelaxed bottom slab has in the high-level database. Defaults to 'bottom_aligned'
    :type bottom_in_name: str, optional

    :param top_out_name: Name the relaxed top slab will have in the high-level database. Defaults to 'top_aligned_relaxed'
    :type top_out_name: str, optional

    :param bottom_out_name: Name the relaxed bottom slab will have in the high-level database. Defaults to 'bottom_aligned_relaxed'
    :type bottom_out_name: str, optional

    :param high_level: Name of the high_level database to use. Defaults to 'True', in which case it is read from the db.json file.
    :type high_level: str or True, optional

    :param prerelax: Whether to run a prerelaxation step before the actual relaxation. Defaults to True.
    :type prerelax: bool, optional

    :param prerelax_calculator: Network potential to use for the prerelaxation. Defaults to 'm3gnet'.
    :type prerelax_calculator: str, optional

    :param prerelax_kwargs: Keyword arguments to pass to the prerelaxation ASE optimizer.
    :type prerelax_kwargs: dict, optional

    :return: FWAction that detours to the relaxation workflow.
    :rtype: FWAction
    """
    _fw_name = "Relax matched slabs"
    required_params = [
        "mp_id_1",
        "mp_id_2",
        "functional",
        "external_pressure",
    ]
    optional_params = [
        "db_file",
        "top_in_name",
        "top_out_name",
        "bottom_in_name",
        "bottom_out_name",
        "high_level",
        "prerelax",
        "prerelax_calculator",
        "prerelax_kwargs",
    ]

    def run_task(self, fw_spec):
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        top_in_name = self.get("top_in_name", "top_aligned")
        top_out_name = self.get("top_out_name", "top_aligned_relaxed")
        bot_in_name = self.get("bottom_in_name", "bottom_aligned")
        bot_out_name = self.get("bottom_out_name", "bottom_aligned_relaxed")
        hl_db = self.get("high_level", True)
        prerelax = self.get("prerelax", True)
        prerelax_calculator = self.get("prerelax_calculator", "m3gnet")
        prerelax_kwargs = self.get("prerelax_kwargs", {})

        name = fw_spec.get("interface_name")
        db = VaspDB(db_file, high_level=hl_db)

        interface_dict = db.find_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
        )

        relaxed_top_present = interface_dict.get(top_out_name)
        relaxed_bot_present = interface_dict.get(bot_out_name)
        comp_params = interface_dict.get("comp_parameters", {})

        inputs = []
        if not relaxed_top_present:
            top_slab = Slab.from_dict(interface_dict[top_in_name])
            top_vis = get_custom_vasp_relax_settings(
                top_slab, comp_params, "slab_pos_relax"
            )
            formula = top_slab.composition.reduced_formula
            miller = "".join(str(s) for s in top_slab.miller_index)
            label_top = (
                "top_slab_"
                + formula
                + miller
                + "_"
                + make_calculation_hash(structure=top_slab, vis=top_vis)
            )
            inputs.append([top_slab, top_vis, label_top])
        if not relaxed_bot_present:
            bot_slab = Slab.from_dict(interface_dict[bot_in_name])
            bot_vis = get_custom_vasp_relax_settings(
                bot_slab, comp_params, "slab_pos_relax"
            )
            formula = bot_slab.composition.reduced_formula
            miller = "".join(str(s) for s in bot_slab.miller_index)
            label_bot = (
                "bot_slab_"
                + formula
                + miller
                + "_"
                + make_calculation_hash(structure=bot_slab, vis=bot_vis)
            )
            inputs.append([bot_slab, bot_vis, label_bot])

        fw_spec["relaxation_inputs"] = inputs

        if inputs:
            wf = dynamic_relax_swf(
                inputs_list=inputs,
                wf_name="Relaxing the matched slabs",
                prerelax_system=prerelax,
                prerelax_calculator=prerelax_calculator,
                prerelax_kwargs=prerelax_kwargs,
                db_file=db_file,
            )
            return FWAction(detours=wf, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_CalcAdhesion(FiretaskBase):
    """Calculate the adhesion from the total energies of static calculations.

    Find the corresponding total energy calculations from the low_level
    database and compute the adhesion energy using the formula:
        E_abs = (top_energy + bot_energy) - inter_energy
    The result is saved in the high_level database in J/m^2.

    :param interface_name: Flag to find the interface entry in the high_level database. Created by the 'interface_name' function in triboflow.firetasks.init_db.pyl
    :type interface_name: str

    :param functional: Functional with which the workflow is run. PBE or SCAN.
    :type functional: str

    :param external_pressure: External pressure in GPa.
    :type external_pressure: float

    :param top_label: Task_label in the low level database to find the total energy calculation of the top slab.
    :type top_label: str

    :param bottom_label: Task_label in the low level database to find the total energy calculation of the bottom slab.
    :type bottom_label: str

    :param interface_label: Task_label in the low level database to find the total energy calculation of the interface.
    :type interface_label: str

    :param db_file: Full path of the db.json file to be used. The default is to use env_chk to find the file.
    :type db_file: str, optional

    :param out_name: Key in the interface collection of the high_level database that contains the calculated adhesion energy as value. Defaults to 'adhesion_energy@min'.
    :type out_name: str, optional

    :param high_level: Name of the high_level database to use. Defaults to 'True', in which case it is read from the db.json file.
    :type high_level: str or True, optional

    :return: FWAction that updates the spec.
    :rtype: FWAction
    """
    _fw_name = "Calculate adhesion energy"
    required_params = [
        "interface_name",
        "functional",
        "top_label",
        "bottom_label",
        "interface_label",
        "external_pressure",
    ]
    optional_params = ["db_file", "out_name", "high_level"]

    def run_task(self, fw_spec):
        name = self.get("interface_name")
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")
        top_label = self.get("top_label")
        bot_label = self.get("bottom_label")
        inter_label = self.get("interface_label")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        out_name = self.get("out_name", "adhesion_energy@min")
        hl_db = self.get("high_level", True)

        db = VaspDB(db_file=db_file)

        top_calc = db.find_data(
            collection="tasks", fltr={"task_label": top_label}
        )
        top_energy = top_calc["output"]["energy"]

        bot_calc = db.find_data(
            collection="tasks", fltr={"task_label": bot_label}
        )
        bot_energy = bot_calc["output"]["energy"]

        inter_calc = db.find_data(
            collection="tasks", fltr={"task_label": inter_label}
        )
        inter_energy = inter_calc["output"]["energy"]
        struct = Structure.from_dict(inter_calc["output"]["structure"])

        area = np.linalg.norm(
            np.cross(struct.lattice.matrix[0], struct.lattice.matrix[1])
        )

        E_abs = (top_energy + bot_energy) - inter_energy

        # Convert adhesion energy from eV/Angstrom^2 to J/m^2
        E_Jm2 = 16.02176565 * E_abs / area

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        db_high.update_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
            new_values={"$set": {out_name: E_Jm2}},
        )

        return FWAction(update_spec=fw_spec)
