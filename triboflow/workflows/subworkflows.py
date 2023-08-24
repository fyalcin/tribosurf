"""SubWorkflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from uuid import uuid4

import numpy as np
from atomate.vasp.fireworks import StaticFW
from atomate.vasp.powerups import add_modify_incar
from fireworks import Workflow, Firework

from hitmen_utils.vasp_tools import (
    get_emin_and_emax,
    get_custom_vasp_static_settings,
)
from hitmen_utils.db_tools import VaspDB
from triboflow.firetasks.PPES import FT_DoPPESCalcs, FT_FitPPES
from triboflow.firetasks.adhesion import FT_CalcAdhesion
from triboflow.firetasks.charge_density_analysis import (
    FT_MakeChargeDensityDiff,
)
from triboflow.firetasks.convergence import FT_Convo
from triboflow.firetasks.dielectric import FT_GetEpsilon
from triboflow.firetasks.structure_manipulation import (
    FT_MakeSlabInDB,
    FT_StartSlabRelax,
    FT_GetRelaxedSlab,
)
from triboflow.firetasks.utils import FT_UpdateCompParams
from triboflow.fireworks.common import run_pes_calc_fw, make_pes_fw
from triboflow.utils.mp_connection import MPConnection


def dielectric_constant_swf(
    structure,
    mpid,
    flag,
    comp_parameters=None,
    spec=None,
    functional="PBE",
    db_file="auto",
    hl_db=True,
    update_bulk=True,
    update_slabs=False,
):
    """
    Subworkflow that calculates dielectric properties and updates the comp_parameters.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which the dielectric properties are calculated.
    mpid : str
        Material Project's material identifier ID of the structure passed.
    flag : str
        Will be used to name the FW of the vasp calc and later to find it in the
        'tasks' collection.
    comp_parameters : dict, optional
        Dictionary of computational parameters. Used to set up the vasp input
        set. The default is {}.
    spec : dict, optional
        fw_spec that can be passed to the FWs. The default is {}.
    functional : str, optional
        Functional for the calculation. Usually SCAN or PBE. Used to select
        the output collection in the high level db. The default is 'PBE'.
    db_file : str, optional
        Path to a db.json file. If 'auto', the standard config folder is used.
        The default is 'auto'.
    hl_db : str or bool, optional
        If a string is given, the high-level database will be chosen based on
        that string. If True, the db.json file will be used to determine the
        name of the high_level database. The default is True.
    update_bulk : bool, optional
        If the bulk entry for the given mpid should be updated.
        The default is True.
    update_slabs : bool, optional
        If the slab entries matching a given mpid should be updated (all miller
        indices). The default is False.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        The dielectric subworkflow.

    """

    # check if epsilon is already calculated for that material.
    if spec is None:
        spec = {}
    if comp_parameters is None:
        comp_parameters = {}
    try:
        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        bulk_data = db_high.find_data(
            collection=f"{functional}.bulk_data", fltr={"mpid": mpid}
        )
        epsilon = bulk_data["comp_parameters"].get("epsilon", False)
    except KeyError:
        epsilon = False

    # if we know that the material is metallic, no need to calculate epsilon
    # since it should be infinite.
    if comp_parameters.get("is_metal", False):
        epsilon = 1000000

    formula = structure.composition.reduced_formula
    wf_name = f"Dielectric calculation WF for {formula} {mpid}"

    if not epsilon:
        vis = get_custom_vasp_static_settings(
            structure, comp_parameters, "bulk_epsilon_from_scratch"
        )
        Calc_Eps_FW = StaticFW(
            structure=structure, name=flag, vasp_input_set=vis
        )

        Get_Eps_FT = FT_GetEpsilon(label=flag, db_file=db_file)

        Update_Data_FT = FT_UpdateCompParams(
            mpid=mpid,
            functional=functional,
            new_params=["epsilon"],
            update_bulk=update_bulk,
            update_slabs=update_slabs,
            db_file=db_file,
            high_level=hl_db,
        )
        Update_FW = Firework(
            tasks=[Get_Eps_FT, Update_Data_FT],
            spec=spec,
            name=flag + "_update_high_level",
        )

        WF = Workflow(
            [Calc_Eps_FW, Update_FW],
            {Calc_Eps_FW.fw_id: [Update_FW.fw_id]},
            name=wf_name,
        )
    else:
        spec.update({"epsilon": epsilon})
        Update_Data_FT = FT_UpdateCompParams(
            mpid=mpid,
            functional=functional,
            new_params=["epsilon"],
            update_bulk=update_bulk,
            update_slabs=update_slabs,
            db_file=db_file,
            high_level=hl_db,
        )
        Update_FW = Firework(
            tasks=[Update_Data_FT], spec=spec, name=flag + "_update_high_level"
        )
        WF = Workflow([Update_FW], name=wf_name)
    return WF


def charge_analysis_swf(
    interface,
    interface_name=None,
    functional="PBE",
    external_pressure=0,
    db_file=None,
    high_level="auto",
    comp_parameters=None,
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
    high_level : str, optional
        Name of the high_level database to use. Defaults to 'auto', in which
        case it is read from the db.json file.
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow to compute the charge density difference of a certain
        interface.

    """
    if comp_parameters is None:
        comp_parameters = {}
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
            high_level=high_level,
        ),
        name=f"Calculate charge density redistribution for {interface_name}",
    )

    SWF = Workflow(
        fireworks=[FW_top, FW_bot, FW_interface, FW_charge_analysis],
        links_dict={
            FW_top.fw_id: [FW_charge_analysis.fw_id],
            FW_bot.fw_id: [FW_charge_analysis.fw_id],
            FW_interface.fw_id: [FW_charge_analysis.fw_id],
        },
        name="Calculate adhesion SWF for {}".format(interface_name),
    )

    return add_modify_incar(SWF)


def adhesion_energy_swf(
    top_slab,
    bottom_slab,
    interface,
    external_pressure=0,
    interface_name=None,
    functional="PBE",
    comp_parameters=None,
    db_file="auto",
    high_level=True,
):
    """Create a subworkflow to compute the adhesion energy for an interface.

    This workflow takes two matched slabs (their cells must be identical) and
    a relaxed interface structure of those slabs and computes the adhesion
    energy. The two matched slabs must be relaxed as well to get correct
    results.
    Output are saved in a high-level database, but may also be written
    as files and copied to a chosen location. Note that this copy operation
    is generally dependent on which machine the calculations are executed,
    and not on the machine where the workflow is submitted. Also, ssh-keys need
    to be set up for remote_copy to work!

    Parameters
    ----------
    db_file
    high_level
    top_slab : pymatgen.core.surface.Slab
        Relaxed top slab of the interface.
    bottom_slab : pymatgen.core.surface.Slab
        Relaxed bottom slab of the interface.
    interface : pymatgen.core.surface.Slab
        Relaxed interface structure.
    external_pressure : float, optional
        External pressure to be applied to the interface. The default is 0.
    interface_name : str, optional
        Unique name to find the interface in the database with.
        The default is None, which will lead to an automatic interface_name
        generation which will be printed on screen.
    functional : str, optional
        Which functional to use; has to be 'PBE' or 'SCAN'. The default is 'PBE'
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow intended to compute the adhesion of a certain interface.

    """
    if comp_parameters is None:
        comp_parameters = {}
    try:
        top_miller = list(top_slab.miller_index)
    except:
        raise AssertionError(
            "You have used {} as an input for <top_slab>.\n"
            "Please use <class 'pymatgen.core.surface.Slab'>"
            " instead.".format(type(top_slab))
        )

    try:
        bot_miller = list(bottom_slab.miller_index)
    except:
        raise AssertionError(
            "You have used {} as an input for <bot_slab>.\n"
            "Please use <class 'pymatgen.core.surface.Slab'>"
            " instead.".format(type(bottom_slab))
        )

    if not interface_name:
        mt = "".join(str(s) for s in top_miller)
        mb = "".join(str(s) for s in bot_miller)
        interface_name = (
            top_slab.composition.reduced_formula
            + "_"
            + mt
            + "_"
            + bottom_slab.composition.reduced_formula
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
        bottom_slab, comp_parameters, "slab_from_scratch"
    )
    vis_interface = get_custom_vasp_static_settings(
        interface, comp_parameters, "slab_from_scratch"
    )

    FW_top = StaticFW(
        structure=top_slab, vasp_input_set=vis_top, name=tag + "top"
    )
    FW_bot = StaticFW(
        structure=bottom_slab, vasp_input_set=vis_bot, name=tag + "bottom"
    )
    FW_interface = StaticFW(
        structure=interface,
        vasp_input_set=vis_interface,
        name=tag + "interface",
    )

    FW_results = Firework(
        FT_CalcAdhesion(
            interface_name=interface_name,
            functional=functional,
            external_pressure=external_pressure,
            top_label=tag + "top",
            bottom_label=tag + "bottom",
            interface_label=tag + "interface",
            db_file=db_file,
            high_level=high_level,
        ),
        name=f"Calculate Adhesion FW for {interface_name}",
    )
    SWF = Workflow(
        fireworks=[FW_top, FW_bot, FW_interface, FW_results],
        links_dict={
            FW_top.fw_id: [FW_results.fw_id],
            FW_bot.fw_id: [FW_results.fw_id],
            FW_interface.fw_id: [FW_results.fw_id],
        },
        name="Calculate adhesion SWF for {}".format(interface_name),
    )

    return add_modify_incar(SWF)


def calc_pes_swf(
    interface,
    interface_name=None,
    functional="PBE",
    pressure=0,
    comp_parameters=None,
    file_output=False,
    output_dir=None,
    remote_copy=False,
    server=None,
    user=None,
    port=None,
    prerelax=True,
    prerelax_calculator="m3gnet",
    prerelax_kwargs=None,
    db_file="auto",
    high_level=True,
):
    """Create a subworkflow to compute the PES for an interface of two slabs.

    This workflow takes two matched slabs (their cells must be identical) as
    input and computes the potential energy surface (PES) for the interface.
    Output are saved in a high-level database, but may also be written
    as files and copied to a chosen location. Note that this copy operation
    is generally dependent on which machine the calculations are executed,
    and not on the machine where the workflow is submitted. Also, ssh-keys need
    to be set up for remote_copy to work!

    Parameters
    ----------
    db_file
    high_level
    interface : pymatgen.core.interface.Interface
        Interface, Slab or Structure object containing the unrelaxed interface.
    interface_name : str, optional
        Unique name to find the interface in the database with.
        The default is None, which will lead to an automatic interface_name
        generation which will be printed on screen.
    functional : str, optional
        Which functional to use; has to be 'PBE' or 'SCAN'. The default is 'PBE'
    pressure : float, optional
        External pressure in GPa to be applied to the interface.
        The default is 0.
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The username on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.
    prerelax : bool, optional
        Whether to perform a prerelaxation using a network potential before starting
        a DFT relaxation. Defaults to True.
    prerelax_calculator : str, optional
        Which network potential to use for the prerelaxation. Defaults to 'm3gnet'.
    prerelax_kwargs : dict, optional
        Keyword arguments to be passed to the ASE calculator for the prerelaxation.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow intended to compute the PES of a certain interface.

    """
    if comp_parameters is None:
        comp_parameters = {}
    if prerelax_kwargs is None:
        prerelax_kwargs = {}
    try:
        top_miller = interface.interface_properties["film_miller"]
        bot_miller = interface.interface_properties["substrate_miller"]
        if not interface_name:
            mt = "".join(str(s) for s in top_miller)
            mb = "".join(str(s) for s in bot_miller)
            interface_name = (
                interface.film.composition.reduced_formula
                + "_"
                + mt
                + "_"
                + interface.substrate.composition.reduced_formula
                + "_"
                + mb
                + "_AutoGen"
            )
            print(
                "\nYour interface name has been automatically generated to be:"
                "\n {}".format(interface_name)
            )
    except:
        if not interface_name:
            mt = mb = "unknown_miller"
            interface_name = (
                interface.film.composition.reduced_formula
                + "_"
                + mt
                + "_"
                + interface.substrate.composition.reduced_formula
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
            "   kpoints_density = 8.5\n"
            "We recommend to pass a comp_parameters dictionary"
            " of the form:\n"
            '   {"use_vdw": <True/False>,\n'
            '    "use_spin": <True/False>,\n'
            '    "is_metal": <True/False>,\n'
            '    "encut": <float>,\n'
            '    "k_dens": <int>}\n'
        )

    tag = interface_name + "_" + str(uuid4())

    db_high = VaspDB(db_file=db_file, high_level=high_level)
    try:
        interface_dict = db_high.find_data(
            collection=f"{functional}.interface_data",
            fltr={"name": interface_name, "pressure": pressure},
        )
        relaxed_structure = interface_dict["relaxed_structure"]
    except KeyError:
        db_high.insert_data(
            collection=functional + ".interface_data",
            data={
                "name": interface_name,
                "pressure": pressure,
                "comp_parameters": comp_parameters,
                "unrelaxed_structure": interface.as_dict(),
                "top_aligned": interface.film.as_dict(),
                "bottom_aligned": interface.substrate.as_dict(),
            },
        )

    FW_1 = run_pes_calc_fw(
        interface=interface,
        interface_name=interface_name,
        external_pressure=pressure,
        functional=functional,
        comp_parameters=comp_parameters,
        tag=tag,
        FW_name="Start PES calcs for " + interface_name,
        prerelax=prerelax,
        prerelax_calculator=prerelax_calculator,
        prerelax_kwargs=prerelax_kwargs,
        db_file=db_file,
        high_level=high_level,
    )

    FW_2 = make_pes_fw(
        interface_name=interface_name,
        functional=functional,
        external_pressure=pressure,
        tag=tag,
        FW_name="Parse PES calcs for " + interface_name,
        file_output=file_output,
        output_dir=output_dir,
        remote_copy=remote_copy,
        server=server,
        user=user,
        port=port,
    )

    SWF = Workflow(
        [FW_1, FW_2],
        {FW_1.fw_id: [FW_2.fw_id]},
        name="Calc PES for " + interface_name + " SWF",
    )
    return SWF


def calc_ppes_swf(
    interface_name,
    functional,
    distance_list=(-0.5, -0.25, 0.0, 0.25, 0.5, 2.5, 3.0, 4.0, 5.0, 7.5),
    out_name="PPES@minimum",
    structure_name="minimum_relaxed",
    spec=None,
    db_file="auto",
    high_level=True,
):
    """
    Generate a subworkflow that calculates a PPES using static calculations.

    For a given interface in the high-level database this subworkflow performs
    static calculations for different distances of the two slabs modeling
    brittle cleavage under mode 1 loading using the rigid separation model.
    The results are saved in an energy vs distance array and saved in the high-level
    database alongside a fit to a UBER curve.

    Parameters
    ----------
    db_file
    high_level
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    distance_list : list of float, optional
        Modification of the equilibrium distance between the slabs.
        The default is [-0.5, -0.25, 0.0, 0.25, 0.5, 2.5, 3.0, 4.0, 5.0, 7.5].
    out_name : str, optional
        Name for the PPES data in the high-level database. The default is
        'PPES@minimum'.
    structure_name : str, optional
        Name of the structure in the interface entry to the high-level database
        for which the PPES should be calculated. The default is
        'minimum_relaxed'.
    spec : dict, optional
        fw_spec that can be passed to the SWF and will be passed on. The
        default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        Subworkflow to calculate the PPES for a certain interface.

    """
    if spec is None:
        spec = {}
    tag = interface_name + "_" + str(uuid4())

    FW_1 = Firework(
        FT_DoPPESCalcs(
            interface_name=interface_name,
            functional=functional,
            distance_list=distance_list,
            tag=tag,
            structure_name=structure_name,
            db_file=db_file,
            high_level=high_level,
        ),
        spec=spec,
        name="PPES Calculations for " + interface_name,
    )

    FW_2 = Firework(
        FT_FitPPES(
            interface_name=interface_name,
            functional=functional,
            distance_list=distance_list,
            out_name=out_name,
            tag=tag,
            db_file=db_file,
            high_level=high_level,
        ),
        spec=spec,
        name="PPES Fitting for " + interface_name,
    )

    SWF = Workflow(
        [FW_1, FW_2],
        {FW_1.fw_id: [FW_2.fw_id]},
        name="Calc PPES for " + interface_name + " SWF",
    )
    return SWF


def make_and_relax_slab_swf(
    bulk_structure,
    miller_index,
    flag,
    comp_parameters=None,
    functional="PBE",
    min_thickness=10.0,
    min_vacuum=15.0,
    relax_type="slab_pos_relax",
    slab_struct_name="unrelaxed_slab",
    out_struct_name="relaxed_slab",
    spec=None,
    file_output=False,
    output_dir=None,
    remote_copy=False,
    server=None,
    user=None,
    port=None,
    print_help=True,
):
    """
    Make and relax a slab.

    Parameters
    ----------
    bulk_structure : pymatgen.core.structure.Structure
        Bulk structure that is used to construct the slab out of.
    miller_index : list of int or str
        Miller indices of the slab to make.
    flag : str
        An identifier to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise, use something
        unique which you can find again.
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.
    functional : str, optional
        Which functional to use; has to be 'PBE' or 'SCAN'. The default is 'PBE'.
    min_thickness : float, optional
        Minimal thickness of the unrelaxed slab in Angstrom. The default is 10.0.
    min_vacuum : float, optional
        Minimal thickness of the vacuum layer in Angstrom. The default is 25.0.
    relax_type : str, optional
        Which type of relaxation to run. See get_custom_vasp_relax_settings from
        hitmen_utils.vasp_tools. The default is 'slab_pos_relax'.
    slab_struct_name : str, optional
        Name of the unrelaxed slab in the high-level database.
        The default is 'unrelaxed_slab'.
    out_struct_name : TYPE, optional
        DESCRIPTION. The default is 'relaxed_slab'.
    spec : dict, optional
        DESCRIPTION. The default is {}.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The username on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow intended to make and relax a slab from a bulk structure.

    """
    if spec is None:
        spec = {}
    if comp_parameters is None:
        comp_parameters = {}
    if type(miller_index) == str:
        miller_str = miller_index
        miller = [int(k) for k in list(miller_index)]
    else:
        miller = miller_index
        miller_str = "".join(str(s) for s in miller_index)

    formula = bulk_structure.composition.reduced_formula

    if flag.startswith("mp-") and flag[3:].isdigit():
        mp_conn = MPConnection()
        formula_from_flag = mp_conn.get_property_from_mp(
            mpid=flag, properties=["formula_pretty"]
        )
        formula_from_flag = formula_from_flag["formula_pretty"]
        if not formula_from_flag == formula:
            raise SystemExit(
                "The chemical formula of your structure ({}) "
                "does not match the chemical formula of the flag "
                "(mp-id) you have chosen which corresponds "
                "to {}.\n".format(formula, formula_from_flag)
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
            '    "encut": <int>,\n'
            '    "k_dens": <int>}\n'
        )

    if print_help:
        db = VaspDB()
        db_file = db.db_file
        print(
            "Once you workflow has finished you can access the "
            "results from the database using this code:\n\n"
            "import pprint\n"
            "from hitmen_utils.db_tools import VaspDB\n"
            f"db = VaspDB(db_file='{db_file}')\n"
            f"results = db.find_data(collection='{functional}.slab_data', fltr={'name': '{flag}'})\n"
            "pprint.pprint(results)\n"
        )

    tag = formula + miller_str + "_" + str(uuid4())

    FTs = []

    FTs.append(
        FT_MakeSlabInDB(
            bulk_structure=bulk_structure,
            miller=miller,
            flag=flag,
            functional=functional,
            min_thickness=min_thickness,
            min_vacuum=min_vacuum,
        )
    )

    FTs.append(
        FT_StartSlabRelax(
            flag=flag,
            miller=miller,
            functional=functional,
            tag=tag,
            comp_parameters=comp_parameters,
            slab_struct_name=slab_struct_name,
            relax_type=relax_type,
        )
    )

    FW = Firework(
        FTs, spec=spec, name="Make and relax " + formula + miller_str + " slab"
    )

    FW2 = Firework(
        FT_GetRelaxedSlab(
            flag=flag,
            miller=miller,
            functional=functional,
            tag=tag,
            struct_out_name=out_struct_name,
            file_output=file_output,
            output_dir=output_dir,
            remote_copy=remote_copy,
            server=server,
            user=user,
            port=port,
        ),
        spec=spec,
        name="Put relaxed " + formula + miller_str + " slab in DB",
    )

    SWF = Workflow(
        [FW, FW2],
        {FW.fw_id: [FW2.fw_id]},
        name="Make and relax " + formula + miller_str + " SWF",
    )
    return SWF


def converge_swf(
    structure,
    conv_type,
    flag,
    comp_parameters=None,
    spec=None,
    functional="PBE",
    deformations=None,
    encut_start=None,
    encut_incr=25,
    k_dens_start=1.0,
    k_dens_incr=0.1,
    k_dens_default=9.0,
    n_converge=3,
    db_file=None,
    file_output=False,
    output_dir=None,
    remote_copy=False,
    server=None,
    user=None,
    port=None,
    print_help=True,
):
    """Subworkflows that converges Encut or kpoints density using fits to an EOS.

    Takes a given structure, computational parameters, and an optional list
    of deformations and uses these deformations to compute a series of
    Birch-Murnaghan equation of state for higher and higher energy cutoffs or
    kpoints density.
    Once bulk modulus and equilibrium volume do not change any longer,
    convergence is reached. Output is printed to the screen and saved in the
    high-level triboflow database where it can be queried using the mp_id
    of the material. Optionally there is also output to a file.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the energy cutoff parameter.
    conv_type : str
        Either "kpoints" or "encut", depending on what to converge.
    flag : str
        An identifier to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise, use something
        unique which you can find again.
    comp_parameters : dict, optional
        Dictionary of computational parameters for the VASP calculations. The
        default is {}.
    spec : dict, optional
        Previous fw_spec that will be updated and/or passed on for child
        Fireworks. The default is {}.
    deformations: list of lists, optional
        List of deformation matrices for the fit to the EOS. Defaults to None,
        which results in 5 volumes from 90% to 110% of the initial volume.
    encut_start : float, optional
        Starting encut value for the first run. Defaults to the largest EMIN
        in the POTCAR.
    encut_incr : float, optional
        Increment for the encut during the convergence. Defaults to 25.
    k_dens_start : float, optional
        Starting kpoint density in 1/Angstrom. Defaults to 1.0
    k_dens_incr : float, optional
        Increment for the kpoint convergence. Can be set quite small since
        there is a check in place to see if a new mesh is actually constructed
        for each density. Defaults to 0.1.
    k_dens_default : float, optional
        Default (quite high) kpoints density for encut convergence studies if
        no k_dens parameter is found in the comp_parameters. The default is 9.0
    n_converge : int, optional
        Number of calculations that have to be inside the convergence
        threshold for convergence to be reached. Defaults to 3.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        None, in which case env_chk will be used in the FT.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The username on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.
    print_help : bool, optional
        Prints a few lines of code that shows how to retrieve the results from
        the database. The default is True.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A subworkflow intended to find the converged ENCUT for a given
        structure.

    """
    if spec is None:
        spec = {}
    if comp_parameters is None:
        comp_parameters = {}
    if conv_type not in ["kpoints", "encut"]:
        raise ValueError(
            '"type" input must be either "kpoints" or'
            '"encut".\nYou have passed {}'.format(conv_type)
        )
    if conv_type == "encut":
        name = (
            "Encut Convergence SWF of " + structure.composition.reduced_formula
        )
        if not encut_start:
            # Get the largest EMIN value of the potcar and round up to the
            # next whole 25.
            vis = get_custom_vasp_static_settings(
                structure, comp_parameters, "bulk_from_scratch"
            )
            encut_dict = get_emin_and_emax(vis.potcar)
            enmax = encut_dict["ENMAX"]
            encut_start = int(25 * np.ceil(enmax / 25))
    elif conv_type == "kpoints":
        name = (
            "Kpoint Convergence SWF of "
            + structure.composition.reduced_formula
        )
    else:
        raise ValueError(
            f'"type" input must be either "kpoints" or'
            '"encut".\nYou have passed {conv_type}.'
        )
    tag = "BM group: {}".format(str(uuid4()))

    if flag.startswith("mp-") and flag[3:].isdigit():
        formula_from_struct = structure.composition.reduced_formula
        mp_conn = MPConnection()
        formula_from_flag = mp_conn.get_property_from_mp(
            mpid=flag, properties=["formula_pretty"]
        )
        formula_from_flag = formula_from_flag["formula_pretty"]

        if not formula_from_flag == formula_from_struct:
            raise SystemExit(
                "The chemical formula of your structure ({}) "
                "does not match the chemical formula of the flag "
                "(mp-id) you have chosen which corresponds "
                "to {}.\n".format(formula_from_struct, formula_from_flag)
            )

    if comp_parameters == {}:
        if conv_type == "encut":
            print(
                "\nNo computational parameters have been defined!\n"
                "Workflow will run with:\n"
                "   ISPIN = 1\n"
                "   ISMEAR = 0\n"
                "   kpoint density kappa = 12.5\n"
                "We recommend to pass a comp_parameters dictionary"
                " of the form:\n"
                '   {"use_vdw": <True/False>,\n'
                '    "use_spin": <True/False>,\n'
                '    "is_metal": <True/False>,\n'
                '    "k_dens": <int>}\n'
            )
        else:
            print(
                "\nNo computational parameters have been defined!\n"
                "Workflow will run with:\n"
                "   ISPIN = 1\n"
                "   ISMEAR = 0\n"
                "   ENCUT = 650\n"
                "We recommend to pass a comp_parameters dictionary"
                " of the form:\n"
                '   {"use_vdw": <True/False>,\n'
                '    "use_spin": <True/False>,\n'
                '    "is_metal": <True/False>,\n'
                '    "encut": <int>}\n'
            )

    if print_help:
        db = VaspDB()
        db_file = db.db_file
        print(
            "Once you workflow has finished you can access the "
            "results from the database using this code:\n\n"
            "import pprint\n"
            "from hitmen_utils.db_tools import VaspDB\n"
            f"db = VaspDB(db_file='{db_file}')\n"
            f"results = db.find_data(collection='{functional}.bulk_data', fltr={'name': '{flag}'})\n"
            "pprint.pprint(results)\n"
        )

    if not comp_parameters.get("functional"):
        comp_parameters["functional"] = functional
    else:
        if not comp_parameters.get("functional") == functional:
            print(
                "The functional set in your computational parameters ({}) "
                "does not match the one given in the input ({})!\n"
                "The functional in the computational parameter has been "
                "overwritten to {}!\n".format(
                    comp_parameters.get("functional"), functional, functional
                )
            )
            comp_parameters["functional"] = functional

    FT_EncutConvo = FT_Convo(
        structure=structure,
        conv_type=conv_type,
        comp_params=comp_parameters,
        tag=tag,
        flag=flag,
        functional=functional,
        deformations=deformations,
        db_file=db_file,
        encut_incr=encut_incr,
        encut_start=encut_start,
        k_dens_start=k_dens_start,
        k_dens_incr=k_dens_incr,
        k_dens_default=k_dens_default,
        n_converge=n_converge,
        file_output=file_output,
        output_dir=output_dir,
        remote_copy=remote_copy,
        server=server,
        user=user,
        port=port,
    )

    FW_C = Firework(FT_EncutConvo, spec=spec, name=name)
    WF = Workflow([FW_C], name=name)

    return WF
