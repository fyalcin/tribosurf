#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 15:36:52 2021

@author: wolloch
"""
import numpy as np
from uuid import uuid4

from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator
from triboflow.utils.structure_manipulation import (
    interface_name,
    slab_from_structure,
)
from hitmen_utils.workflows import dynamic_relax_swf
from hitmen_utils.vasp_tools import get_custom_vasp_relax_settings


@explicit_serialize
class FT_RetrieveMatchedSlabs(FiretaskBase):
    """Retrieve relaxed matched slabs and save them in the high_level database.

    Get the relaxed aligned top and bottom slabs from the low level database
    and save them in the interface_data collection of the high level database.

    Parameters
    ----------
    mp_id_1 : str
        MaterialsProject ID number for the first material
    mp_id_2 : str
        MaterialsProject ID number for the second material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    external_pressure : float
        External pressure in GPa.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    top_out_name : str, optional
        Name the relaxed top slab will have in the high-level database.
        Defaults to 'top_aligned_relaxed'
    bottom_out_name : str, optional
        Name the relaxed bottom slab will have in the high-level database.
        Defaults to 'bottom_aligned_relaxed'
    high_level : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    """

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
        pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        top_out_name = self.get("top_out_name", "top_aligned_relaxed")
        bot_out_name = self.get("bottom_out_name", "bottom_aligned_relaxed")
        hl_db = self.get("high_level", True)

        name = fw_spec.get("interface_name")
        input_list = fw_spec.get("relaxation_inputs")

        if input_list:
            nav = Navigator(db_file)
            for i in input_list:
                label = i[-1]
                miller = i[0].miller_index
                calc = nav.find_data(
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
                nav_high = Navigator(db_file, high_level=hl_db)
                nav_high.update_data(
                    collection=functional + ".interface_data",
                    fltr={"name": name, "pressure": pressure},
                    new_values={"$set": {out_name: slab.as_dict()}},
                )

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_RelaxMatchedSlabs(FiretaskBase):
    """Start the relaxation of the matched slabs.

    Get the aligned top and bottom slabs from the high level database and
    start their relaxation runs.

    Parameters
    ----------
    mp_id_1 : str
        MaterialsProject ID number for the first material
    mp_id_2 : str
        MaterialsProject ID number for the second material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    external_pressure : float
        External pressure in GPa.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    top_in_name : str, optional
        Name the unrelaxed top slab has in the high-level database.
        Defaults to 'top_aligned'
    bottom_in_name : str, optional
        Name the unrelaxed bottom slab has in the high-level database.
        Defaults to 'bottom_aligned'
    top_out_name : str, optional
        Name the relaxed top slab will have in the high-level database.
        Defaults to 'top_aligned_relaxed'
    bottom_out_name : str, optional
        Name the relaxed bottom slab will have in the high-level database.
        Defaults to 'bottom_aligned_relaxed'
    high_level : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    prerelax : bool, optional
        Whether to run a prerelaxation step before the actual relaxation.
        Defaults to True.
    prerelax_calculator : str, optional
        Network potential to use for the prerelaxation. Defaults to 'm3gnet'.
    prerelax_kwargs : dict, optional
        Keyword arguments to pass to the prerelaxation ASE optimizer.
    """

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
        pressure = self.get("external_pressure")
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
        nav = Navigator(db_file=db_file, high_level=hl_db)

        interface_dict = nav.find_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "pressure": pressure},
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
            label = "top_slab_" + formula + miller + "_" + str(uuid4())
            inputs.append([top_slab, top_vis, label])
        if not relaxed_bot_present:
            bot_slab = Slab.from_dict(interface_dict[bot_in_name])
            bot_vis = get_custom_vasp_relax_settings(
                bot_slab, comp_params, "slab_pos_relax"
            )
            formula = bot_slab.composition.reduced_formula
            miller = "".join(str(s) for s in bot_slab.miller_index)
            label = "bot_slab_" + formula + miller + "_" + str(uuid4())
            inputs.append([bot_slab, bot_vis, label])

        fw_spec["relaxation_inputs"] = inputs

        if inputs:
            wf = dynamic_relax_swf(
                inputs_list=inputs,
                wf_name="Relaxing the matched slabs",
                prerelax_system=prerelax,
                prerelax_calculator=prerelax_calculator,
                prerelax_kwargs=prerelax_kwargs,
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

    Parameters
    ----------
    interface_name : str
        Flag to find the interface entry in the high_level database. Created
        by the 'interface_name' function in triboflow.firetasks.init_db.pyl
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    external_pressure : float
        External pressure in GPa.
    top_label : str
        Task_label in the low level database to find the total energy
        calculation of the top slab.
    bottom_label : str
        Task_label in the low level database to find the total energy
        calculation of the bottom slab.
    interface_label : str
        Task_label in the low level database to find the total energy
        calculation of the interface.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    out_name: str, optional
        Key in the interface collection of the high_level database that
        contains the calculated adhesion energy as value. Defaults to
        'adhesion_energy@min'.
    high_level : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    """

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
        pressure = self.get("external_pressure")
        top_label = self.get("top_label")
        bot_label = self.get("bottom_label")
        inter_label = self.get("interface_label")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        out_name = self.get("out_name", "adhesion_energy@min")
        hl_db = self.get("high_level", True)

        nav = Navigator(db_file=db_file)

        top_calc = nav.find_data(
            collection="tasks", fltr={"task_label": top_label}
        )
        top_energy = top_calc["output"]["energy"]

        bot_calc = nav.find_data(
            collection="tasks", fltr={"task_label": bot_label}
        )
        bot_energy = bot_calc["output"]["energy"]

        inter_calc = nav.find_data(
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

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "pressure": pressure},
            new_values={"$set": {out_name: E_Jm2}},
        )

        return FWAction(update_spec=fw_spec)
