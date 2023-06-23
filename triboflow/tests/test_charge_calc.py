#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 16:23:42 2021

@author: mwo
"""
from uuid import uuid4

from fireworks import LaunchPad
from fireworks import Workflow, Firework

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.fireworks.core import StaticFW

from triboflow.utils.vasp_tools import get_custom_vasp_static_settings
from triboflow.utils.database import NavigatorMP
from triboflow.phys.interface_matcher import InterfaceMatcher
from triboflow.phys.shaper import Shaper
from triboflow.firetasks.charge_density_analysis import (
    FT_MakeChargeDensityDiff,
)


def charge_analysis_swf(
    interface,
    interface_name=None,
    functional="PBE",
    external_pressure=0,
    db_file=None,
    high_level_db="auto",
    comp_parameters={},
):
    """Subworkflow to compute the charge redistribution of an interface.

    This workflow takes an interface object as an input and makes 3 static
    calculations of the interface, and the top and bottom slabs, respectively.
    The charge densities are parsed and subsequently analyzed. The difference
    is computed (chgcar_interface - chgcar_top - chgcar_bottom) and charge
    redistribution is computed similar to PRL 121, 026804 (2018).

    Output are saved in a high-level database.

    Parameters
    ----------.
    interface : pymatgen.core.interface.Interface
        Relaxed interface structure.
    interface_name : str, optional
        Unique name to find the interface in the database with.
        The default is None, which will lead to an automatic interface_name
        generation which will be printed on screen.
    functional : str, optional
        Which functional to use; has to be 'PBE' or 'SCAN'. The default is 'PBE'
    external_pressure : float, optional
        External pressure to be applied to the interface. The default is 0.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    high_level_db : str, optional
        Name of the high_level database to use. Defaults to 'auto', in which
        case it is read from the db.json file.
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow to compute the charge density differnece of a certain
        interface.

    """
    top_slab = interface.film
    bot_slab = interface.substrate
    try:
        top_miller = interface.interface_properties["film_miller"]
        bot_miller = interface.interface_properties["substrate_miller"]
    except:
        top_miller = bot_miller = [0, 0, 0]
        RuntimeWarning(
            'Interface object is missing the "interface_properties"\n'
            '"film_miller", and "substrate_miller"'
        )

    if not interface_name:
        mt = "".join(str(s) for s in top_miller)
        mb = "".join(str(s) for s in bot_miller)
        interface_name = (
            top_slab.composition.reduced_formula
            + "_"
            + mt
            + "_"
            + bot_slab.composition.reduced_formula
            + "_"
            + mb
            + "_AutoGen"
        )
        print(
            "\nYour interface name has been automatically generated to be:"
            "\n {}".format(interface_name)
        )

    if comp_parameters == {}:
        print(
            "\nNo computational parameters have been defined!\n"
            "Workflow will run with:\n"
            "   ISPIN = 1\n"
            "   ISMEAR = 0\n"
            "   ENCUT = 520\n"
            "   kpoint density kappa = 5000\n"
            "We recommend to pass a comp_parameters dictionary"
            " of the form:\n"
            '   {"use_vdw": <True/False>,\n'
            '    "use_spin": <True/False>,\n'
            '    "is_metal": <True/False>,\n'
            '    "encut": <float>,\n'
            '    "k_dens": <int>}\n'
        )

    tag = interface_name + "_" + str(uuid4())

    vis_top = get_custom_vasp_static_settings(
        top_slab, comp_parameters, "slab_from_scratch"
    )
    vis_bot = get_custom_vasp_static_settings(
        bot_slab, comp_parameters, "slab_from_scratch"
    )
    vis_interface = get_custom_vasp_static_settings(
        interface, comp_parameters, "slab_from_scratch"
    )

    FW_top = StaticFW(
        structure=top_slab,
        vasp_input_set=vis_top,
        name=tag + "top",
        vasptodb_kwargs={"store_volumetric_data": ["chgcar"]},
    )
    FW_bot = StaticFW(
        structure=bot_slab,
        vasp_input_set=vis_bot,
        name=tag + "bottom",
        vasptodb_kwargs={"store_volumetric_data": ["chgcar"]},
    )
    FW_interface = StaticFW(
        structure=interface,
        vasp_input_set=vis_interface,
        name=tag + "interface",
        vasptodb_kwargs={"store_volumetric_data": ["chgcar"]},
    )

    FW_charge_analysis = Firework(
        FT_MakeChargeDensityDiff(
            interface=interface,
            interface_name=interface_name,
            interface_calc_name=tag + "interface",
            top_calc_name=tag + "top",
            bot_calc_name=tag + "bottom",
            functional=functional,
            external_pressure=external_pressure,
            db_file=db_file,
            high_level_db=high_level_db,
        ),
        name=f"Calculate charge density redistribution for {interface_name}",
    )

    SWF = Workflow(
        fireworks=[FW_top, FW_bot, FW_interface, FW_charge_analysis],
        links_dict={
            FW_top: [FW_charge_analysis],
            FW_bot: [FW_charge_analysis],
            FW_interface: [FW_charge_analysis],
        },
        name="Calculate adhesion SWF for {}".format(interface_name),
    )

    return add_modify_incar(SWF)


if __name__ == "__main__":
    comp_params = {
        "encut": 400,
        "use_spin": False,
        "k_dens": 3.0,
        "is_metal": True,
        "epsilon": 100000,
        "use_vdw": False,
        "functional": "PBE",
    }
    nav_mp = NavigatorMP()
    graphite, _ = nav_mp.get_low_energy_structure("C", mp_id="mp-48")
    nickel, _ = nav_mp.get_low_energy_structure("Ni", mp_id="mp-23")

    gr_conv = SpacegroupAnalyzer(
        graphite
    ).get_conventional_standard_structure()
    ni_conv = SpacegroupAnalyzer(nickel).get_conventional_standard_structure()

    sg_params_ni = {
        "miller": [1, 1, 1],
        "slab_thick": 4,
        "vac_thick": 20,
        "max_normal_search": 1,
        "lll_reduce": True,
        "primitive": True,
        "tol": 0.1,
    }
    sg_params_gr = {
        "miller": [0, 0, 1],
        "slab_thick": 3,
        "vac_thick": 20,
        "max_normal_search": 1,
        "lll_reduce": True,
        "primitive": True,
        "tol": 0.1,
    }

    gr_slab_d, _ = Shaper.generate_slabs(gr_conv, sg_params_gr)
    ni_slab_d, _ = Shaper.generate_slabs(ni_conv, sg_params_ni)
    gr_slab = gr_slab_d[(0, 0, 1)][0]
    ni_slab = ni_slab_d[(1, 1, 1)][0]

    IM = InterfaceMatcher(gr_slab, ni_slab)
    interface = IM.get_interface()

    pressure = 0.0

    WF = charge_analysis_swf(
        interface=interface,
        comp_parameters=comp_params,
        high_level_db="test",
        functional="PBE",
        external_pressure=pressure,
    )
    lpad = LaunchPad.auto_load()
    lpad.add_wf(WF)
    # FW1 = Firework(tasks=[FT_MakeChargeCalc(structure=slab,
    #                                        comp_params=comp_params,
    #                                        calc_name=calc_name)],
    #               name='compute_charge_FW')
    # FW2 = Firework(tasks=[FT_GetCharge(calc_name=calc_name)],
    #               name='get_charge_FW')

    # WF = Workflow([FW1, FW2], {FW1: [FW2]}, name='Test_WF_Charge')
    # lpad = LaunchPad.auto_load()
    # lpad.add_wf(WF)
