#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 16:22:56 2022

@author: wolloch
"""

import numpy as np
from atomate.utils.utils import env_chk
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.powerups import add_modify_incar
from fireworks import Workflow, FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.interface import Interface
from pymatgen.io.vasp.outputs import Chgcar
from scipy.integrate import romb, simpson

from htflow_utils.db_tools import VaspDB
from htflow_utils.vasp_tools import get_custom_vasp_static_settings


def plot_charge_profile(
    chgcar: Chgcar, axis: int = 2, xmin: float = None, xmax: float = None
) -> None:
    """Plots the charge density profile along the given axis.

    :param chgcar: Charge density object.
    :type chgcar: pymatgen.io.vasp.outputs.Chgcar

    :param axis: Axis along which the profile is plotted.
    :type axis: int

    :param xmin: Minimum value of the x-axis.
    :type xmin: float

    :param xmax: Maximum value of the x-axis.
    :type xmax: float

    :return: None
    :rtype: NoneType
    """
    x = chgcar.get_axis_grid(axis)
    y = chgcar.get_average_along_axis(axis)
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = plt.axes()
    ax.plot(x, y)
    if xmin or xmax:
        delta = 1.5 * (xmax - xmin)
        ax.set_xlim(xmin - delta, xmax + delta)
        plt.axvline(x=xmin, color="black", linestyle=":")
        plt.axvline(x=xmax, color="black", linestyle=":")
    fig.savefig(f"charge_profile_{chgcar.structure.formula}.png")


def get_interface_region(interface: Interface) -> tuple:
    """Returns the z-coordinates of the interface region.

    :param interface: Interface object.
    :type interface: pymatgen.core.interface.Interface

    :return: zmin, zmax
    :rtype: tuple
    """
    z_sub = [
        s.coords[2]
        for s in interface.sites
        if s.properties["interface_label"] == "substrate"
    ]
    z_film = [
        s.coords[2]
        for s in interface.sites
        if s.properties["interface_label"] == "film"
    ]
    return max(z_sub), min(z_film)


def make_charge_differences(
    interface: Interface,
    chgcar_int: Chgcar,
    chgcar_bot: Chgcar,
    chgcar_top: Chgcar,
) -> dict:
    """Computes the charge density difference between the interface and the
    top and bottom slabs.

    :param interface: Interface object.
    :type interface: pymatgen.core.interface.Interface

    :param chgcar_int: Charge density of the interface.
    :type chgcar_int: pymatgen.io.vasp.outputs.Chgcar

    :param chgcar_bot: Charge density of the bottom slab.
    :type chgcar_bot: pymatgen.io.vasp.outputs.Chgcar

    :param chgcar_top: Charge density of the top slab.
    :type chgcar_top: pymatgen.io.vasp.outputs.Chgcar

    :return: Dictionary containing the charge density differences.
    :rtype: dict
    """
    chgcar_diff = chgcar_int - chgcar_top - chgcar_bot
    abs_chgcar_diff = chgcar_diff.copy()
    abs_chgcar_diff.data["total"] = abs(abs_chgcar_diff.data["total"])

    zmin, zmax = get_interface_region(interface)
    z = chgcar_diff.get_axis_grid(2)
    dz = z[1] - z[0]

    profile = chgcar_diff.get_average_along_axis(2)
    profile_inter_region = np.asarray(
        [
            profile[i]
            for i in np.arange(chgcar_diff.dim[2])
            if zmin <= z[i] <= zmax
        ]
    )

    abs_profile = abs_chgcar_diff.get_average_along_axis(2)
    abs_profile_inter_region = np.asarray(
        [
            abs_profile[i]
            for i in np.arange(abs_chgcar_diff.dim[2])
            if zmin <= z[i] <= zmax
        ]
    )

    try:
        rho_total = romb(y=abs(profile), dx=dz) / (zmax - zmin)
        rho_total_abs = romb(y=abs_profile, dx=dz) / (zmax - zmin)
    except:
        rho_total = simpson(y=abs(profile), dx=dz) / (zmax - zmin)
        rho_total_abs = simpson(y=abs_profile, dx=dz) / (zmax - zmin)

    try:
        rho_inter_region = romb(y=abs(profile_inter_region), dx=dz) / (
            zmax - zmin
        )
        rho_inter_region_abs = romb(y=abs_profile_inter_region, dx=dz) / (
            zmax - zmin
        )
    except:
        rho_inter_region = simpson(y=abs(profile_inter_region), dx=dz) / (
            zmax - zmin
        )
        rho_inter_region_abs = simpson(y=abs_profile_inter_region, dx=dz) / (
            zmax - zmin
        )

    return {
        "rho_classic": rho_inter_region,
        "rho_all": rho_total,
        "rho_new": rho_inter_region_abs,
        "rho_all_new": rho_total_abs,
    }


@explicit_serialize
class FT_MakeChargeCalc(FiretaskBase):
    """Firetask to make a charge density calculation.

    :param structure: Structure object.
    :type structure: pymatgen.core.structure.Structure

    :param comp_params: Dictionary containing computational parameters.
    :type comp_params: dict

    :param calc_name: Name of the calculation.
    :type calc_name: str

    :param db_file: Path to the database file.
    :type db_file: str

    :param high_level: Name of the high-level database.
    :type high_level: str

    :return: FWAction that detours to a StaticFW.
    :rtype: FWAction
    """

    _fw_name = "Make charge density calculation"
    required_params = ["structure", "comp_params", "calc_name"]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        struct = self.get("structure")
        comp_params = self.get("comp_params")
        label = self.get("calc_name")

        vis = get_custom_vasp_static_settings(
            struct, comp_params, "slab_from_scratch"
        )

        FW = StaticFW(
            structure=struct,
            vasp_input_set=vis,
            name=label,
            vasptodb_kwargs={"store_volumetric_data": ["chgcar"]},
        )
        WF = add_modify_incar(Workflow.from_Firework(FW))
        return FWAction(detours=WF)


@explicit_serialize
class FT_MakeChargeDensityDiff(FiretaskBase):
    """Firetask to make a charge density difference calculation.

    :param interface: Interface object.
    :type interface: pymatgen.core.interface.Interface

    :param interface_name: Name of the interface.
    :type interface_name: str

    :param interface_calc_name: Name of the interface calculation.
    :type interface_calc_name: str

    :param top_calc_name: Name of the top slab calculation.
    :type top_calc_name: str

    :param bot_calc_name: Name of the bottom slab calculation.
    :type bot_calc_name: str

    :param functional: Functional to be used.
    :type functional: str

    :param external_pressure: External pressure in GPa.
    :type external_pressure: float

    :param db_file: Path to the database file.
    :type db_file: str

    :param high_level: Name of the high-level database.
    :type high_level: str

    :return: None
    :rtype: NoneType
    """

    _fw_name = "Make charge density difference"
    required_params = [
        "interface",
        "interface_name",
        "interface_calc_name",
        "top_calc_name",
        "bot_calc_name",
        "functional",
        "external_pressure",
    ]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        interface = self.get("interface")
        name = self.get("interface_name")
        label_int = self.get("interface_calc_name")
        label_top = self.get("top_calc_name")
        label_bot = self.get("bot_calc_name")
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        hl_db = self.get("high_level", "auto")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        db = VaspDB(db_file=db_file)
        chgcar_int = db.get_chgcar_from_label(label_int)
        chgcar_top = db.get_chgcar_from_label(label_top)
        chgcar_bot = db.get_chgcar_from_label(label_bot)

        rho_dict = make_charge_differences(
            interface, chgcar_int, chgcar_bot, chgcar_top
        )

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        db_high.update_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
            new_values={"$set": {"charge_density_redist": rho_dict}},
            upsert=True,
        )

        return
