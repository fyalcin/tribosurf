#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 16:22:56 2022

@author: wolloch
"""

import numpy as np
from scipy.integrate import romb, simpson

from fireworks import Workflow, FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.fireworks.core import StaticFW
from atomate.utils.utils import env_chk

from hitmen_utils.vasp_tools import get_custom_vasp_static_settings
from triboflow.utils.database import Navigator


def plot_charge_profile(chgcar, axis=2, xmin=None, xmax=None):
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


def get_interface_region(interface):
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


def make_charge_differences(interface, chgcar_int, chgcar_bot, chgcar_top):
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
        rho_inter_region = romb(y=abs(profile_inter_region), dx=dz) / (zmax - zmin)
        rho_inter_region_abs = romb(y=abs_profile_inter_region, dx=dz) / (zmax - zmin)
    except:
        rho_inter_region = simpson(y=abs(profile_inter_region), dx=dz) / (zmax - zmin)
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
    required_params = ["structure", "comp_params", "calc_name"]
    optional_params = ["db_file", "high_level_db"]

    def run_task(self, fw_spec):
        struct = self.get("structure")
        comp_params = self.get("comp_params")
        label = self.get("calc_name")

        vis = get_custom_vasp_static_settings(struct, comp_params, "slab_from_scratch")

        FW = StaticFW(
            structure=struct,
            vasp_input_set=vis,
            name=label,
            vasptodb_kwargs={"store_volumetric_data": ["chgcar"]},
        )
        WF = add_modify_incar(Workflow.from_Firework(FW))
        return FWAction(detours=WF)


@explicit_serialize
class FT_GetCharge(FiretaskBase):
    required_params = ["calc_name"]
    optional_params = ["db_file"]

    def run_task(self, fw_spec):
        label = self.get("calc_name")
        if self.get("db_file", None):
            nav = Navigator(db_file=self.get("db_file"))
        else:
            nav = Navigator()
        chgcar = nav.get_chgcar_from_label(label)
        plot_charge_profile(chgcar)


@explicit_serialize
class FT_MakeChargeDensityDiff(FiretaskBase):
    required_params = [
        "interface",
        "interface_name",
        "interface_calc_name",
        "top_calc_name",
        "bot_calc_name",
        "functional",
        "external_pressure",
    ]
    optional_params = ["db_file", "high_level_db"]

    def run_task(self, fw_spec):
        interface = self.get("interface")
        name = self.get("interface_name")
        label_int = self.get("interface_calc_name")
        label_top = self.get("top_calc_name")
        label_bot = self.get("bot_calc_name")
        functional = self.get("functional")
        pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        hl_db = self.get("high_level_db", "auto")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        nav = Navigator(db_file=db_file)
        chgcar_int = nav.get_chgcar_from_label(label_int)
        chgcar_top = nav.get_chgcar_from_label(label_top)
        chgcar_bot = nav.get_chgcar_from_label(label_bot)

        rho_dict = make_charge_differences(
            interface, chgcar_int, chgcar_bot, chgcar_top
        )

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "pressure": pressure},
            new_values={"$set": {"charge_density_redist": rho_dict}},
            upsert=True,
        )

        return
