"""Main Workflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from fireworks import Workflow, Firework

from triboflow.firetasks.adhesion import (
    FT_RelaxMatchedSlabs,
    FT_RetrieveMatchedSlabs,
)
from triboflow.firetasks.run_slabs_wfs import (
    GetCandidatesForHeteroStructure,
)
from triboflow.firetasks.start_swfs import (
    FT_StartAdhesionSWF,
    FT_StartBulkConvoSWF,
    # FT_StartDielectricSWF,
    FT_StartPESCalcSWF,
    FT_StartChargeAnalysisSWF,
)
from triboflow.firetasks.structure_manipulation import (
    FT_AddBulkToDB,
    FT_MakeHeteroStructure,
    FT_StartBulkPreRelax,
)
from triboflow.utils.mp_connection import material_from_mp
from triboflow.utils.check_WF_inputs import check_hetero_wf_inputs


def optimize_bulk_wf(inputs):
    """Return main workflow for bulk optimization within Triboflow.

    Parameters
    ----------
    inputs : dict
        Dictionary containing sub-dictionaries with material information and
        parameters for the interface creation and computational settings.

    Returns
    -------
    wf_list : FireWorks Workflow
        Main Triboflow workflow for bulk optimization.

    """

    material = inputs["material"]
    computational_params = inputs["computational_params"]
    db_file = inputs["database_params"]["db_file"]
    high_level = inputs["database_params"]["high_level"]

    struct, mpid = material_from_mp(material)

    formula = struct.composition.reduced_formula

    functional = computational_params.get("functional", "PBE")

    wf_list = []

    add_bulk = Firework(
        FT_AddBulkToDB(
            mpid=mpid,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
            custom_data={"comp_parameters": computational_params},
        ),
        name=f"Add {formula} ({mpid}) to DB",
    )
    wf_list.append(add_bulk)

    pre_relaxation = Firework(
        FT_StartBulkPreRelax(
            mp_id=mpid,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
            k_dense=computational_params.get("k_dense", 9)
        ),
        name=f"Start pre-relaxation for {formula} ({mpid})",
        parents=[add_bulk],
    )
    wf_list.append(pre_relaxation)

    converge_encut = Firework(
        FT_StartBulkConvoSWF(
            conv_type="encut",
            mp_id=mpid,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
        ),
        name=f"Start encut convergence for {formula} ({mpid})",
        parents=[pre_relaxation],
    )
    wf_list.append(converge_encut)

    converge_kpoints = Firework(
        FT_StartBulkConvoSWF(
            conv_type="kpoints",
            mp_id=mpid,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
        ),
        name=f"Start kpoints convergence for {formula} ({mpid})",
        parents=[converge_encut],
    )
    wf_list.append(converge_kpoints)

    wf_name = f"Converge_bulk_{formula}_{functional}"

    wf = Workflow(wf_list, name=wf_name)

    return wf


def heterogeneous_wf(inputs):
    """Return main workflow for heterogeneous interfaces within Triboflow.

    Parameters
    ----------
    inputs : dict
        Dictionary containing sub-dictionaries with material information and
        parameters for the interface creation and computational settings.

    Returns
    -------
    fw_list : FireWorks Workflow
        Main Triboflow workflow for heterogeneous interfaces.

    """

    inputs = check_hetero_wf_inputs(inputs)

    material_1 = inputs["material_1"]
    material_2 = inputs["material_2"]
    computational_params = inputs["computational_params"]
    interface_params = inputs["interface_params"]
    sg_params_1 = inputs["sg_params_1"]
    sg_params_2 = inputs["sg_params_2"]
    sg_filter_1 = inputs["sg_filter_1"]
    sg_filter_2 = inputs["sg_filter_2"]
    db_file = inputs["database_params"]["db_file"]
    high_level = inputs["database_params"]["high_level"]

    struct_1, mpid_1 = material_from_mp(material_1)
    if material_1["mpid"] and mpid_1 != material_1["mpid"]:
        raise ValueError(
            f"MPID {material_1['mpid']} does not match the MPID\n"
            f"for the formula {material_1['formula']} returned by the Materials Project API\n"
            "Please update your inputs with matching formulas and mpids."
        )
    elif struct_1.composition.reduced_formula != material_1["formula"]:
        raise ValueError(
            f"Formula '{material_1['formula']}' does not match the formula\n"
            f"for the MPID {material_1['mpid']} returned by the Materials Project API:\n"
            f"'{struct_1.composition.reduced_formula}'. Please update your inputs with matching formulas and mpids."
        )
    struct_2, mpid_2 = material_from_mp(material_2)
    if material_2["mpid"] and mpid_2 != material_2["mpid"]:
        raise ValueError(
            f"MPID {material_2['mpid']} does not match the MPID\n"
            f"for the formula {material_2['formula']} returned by the Materials Project API\n"
            "Please update your inputs with matching formulas and mpids."
        )
    elif struct_2.composition.reduced_formula != material_2["formula"]:
        raise ValueError(
            f"Formula '{material_2['formula']}' does not match the formula\n"
            f"for the MPID {material_2['mpid']} returned by the Materials Project API:\n"
            f"'{struct_2.composition.reduced_formula}'. Please update your inputs with matching formulas and mpids."
        )


    formula_1 = struct_1.composition.reduced_formula
    formula_2 = struct_2.composition.reduced_formula

    functional = computational_params.get("functional", "PBE")

    fw_list = []

    # external_pressure might default to None, so we have to check for that
    external_pressure = (
        interface_params.get("external_pressure", 0.0) or 0.0
    )

    add_bulk_m1 = Firework(
        FT_AddBulkToDB(
            mpid=mpid_1,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
            custom_data={"comp_parameters": computational_params},
        ),
        name=f"Add {formula_1} ({mpid_1}) to DB",
    )
    fw_list.append(add_bulk_m1)

    pre_relaxation_m1 = Firework(
        FT_StartBulkPreRelax(
            mp_id=mpid_1,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
        ),
        name=f"Start pre-relaxation for {formula_1} ({mpid_1})",
        parents=[add_bulk_m1],
    )
    fw_list.append(pre_relaxation_m1)

    converge_encut_m1 = Firework(
        FT_StartBulkConvoSWF(
            conv_type="encut",
            mp_id=mpid_1,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
        ),
        name=f"Start encut convergence for {formula_1} ({mpid_1})",
        parents=[pre_relaxation_m1],
    )
    fw_list.append(converge_encut_m1)

    converge_kpoints_m1 = Firework(
        FT_StartBulkConvoSWF(
            conv_type="kpoints",
            mp_id=mpid_1,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
        ),
        name=f"Start kpoints convergence for {formula_1} ({mpid_1})",
        parents=[converge_encut_m1],
    )
    fw_list.append(converge_kpoints_m1)

    final_params_parents = [converge_kpoints_m1]

    if mpid_2 != mpid_1:
        add_bulk_m2 = Firework(
            FT_AddBulkToDB(
                mpid=mpid_2,
                functional=functional,
                db_file=db_file,
                high_level=high_level,
                custom_data={"comp_parameters": computational_params},
            ),
            name=f"Add {formula_2} ({mpid_2}) to DB",
        )
        fw_list.append(add_bulk_m2)

        pre_relaxation_m2 = Firework(
            FT_StartBulkPreRelax(
                mp_id=mpid_2,
                functional=functional,
                db_file=db_file,
                high_level=high_level,
            ),
            name=f"Start pre-relaxation for {formula_2} ({mpid_2})",
            parents=[add_bulk_m2],
        )
        fw_list.append(pre_relaxation_m2)

        converge_encut_m2 = Firework(
            FT_StartBulkConvoSWF(
                conv_type="encut",
                mp_id=mpid_2,
                functional=functional,
                db_file=db_file,
                high_level=high_level,
            ),
            name=f"Start encut convergence for {formula_2} ({mpid_2})",
            parents=[pre_relaxation_m2],
        )
        fw_list.append(converge_encut_m2)

        converge_kpoints_m2 = Firework(
            FT_StartBulkConvoSWF(
                conv_type="kpoints",
                mp_id=mpid_2,
                functional=functional,
                db_file=db_file,
                high_level=high_level,
            ),
            name=f"Start kpoints convergence for {formula_2} ({mpid_2})",
            parents=[converge_encut_m2],
        )
        fw_list.append(converge_kpoints_m2)

        final_params_parents.append(converge_kpoints_m2)

    bulk_coll = f"{functional}.bulk_data"
    surfen_coll = f"{functional}.surfen_data"

    find_candidate_slabs = Firework(
        GetCandidatesForHeteroStructure(
            mpid_1=mpid_1,
            mpid_2=mpid_2,
            comp_params_1={},
            comp_params_2={},
            interface_params=interface_params,
            sg_params_1=sg_params_1,
            sg_params_2=sg_params_2,
            sg_filter_1=sg_filter_1,
            sg_filter_2=sg_filter_2,
            db_file=db_file,
            high_level=high_level,
            surfen_coll=surfen_coll,
            bulk_coll=bulk_coll,
            add_full_relax=True,
        ),
        name="Get candidate slabs for heterostructure generation",
        parents=final_params_parents,
    )

    fw_list.append(find_candidate_slabs)

    make_interface = Firework(
        FT_MakeHeteroStructure(
            mp_id_1=mpid_1,
            mp_id_2=mpid_2,
            interface_params=interface_params,
            functional=functional,
            db_file=db_file,
            high_level=high_level,
        ),
        name="Match the interface",
        parents=[find_candidate_slabs],
    )
    fw_list.append(make_interface)

    relax_matched_slabs = Firework(
        FT_RelaxMatchedSlabs(
            mp_id_1=mpid_1,
            mp_id_2=mpid_2,
            functional=functional,
            external_pressure=external_pressure,
            prerelax=True,
            db_file=db_file,
            high_level=high_level,
        ),
        name="Fully relax the matched slabs",
        parents=[make_interface],
    )
    fw_list.append(relax_matched_slabs)

    calc_pes_points = Firework(
        FT_StartPESCalcSWF(
            mp_id_1=mpid_1,
            mp_id_2=mpid_2,
            functional=functional,
            external_pressure=external_pressure,
            prerelax=True,
            db_file=db_file,
            high_level=high_level,
        ),
        name="Compute PES high-symmetry points",
        parents=[make_interface],
    )
    fw_list.append(calc_pes_points)

    retrieve_matched_slabs = Firework(
        FT_RetrieveMatchedSlabs(
            mp_id_1=mpid_1,
            mp_id_2=mpid_2,
            functional=functional,
            external_pressure=external_pressure,
            db_file=db_file,
            high_level=high_level,
        ),
        name="Retrieve relaxed matched slabs",
        parents=[relax_matched_slabs],
    )
    fw_list.append(retrieve_matched_slabs)

    compute_adhesion = Firework(
        FT_StartAdhesionSWF(
            mp_id_1=mpid_1,
            mp_id_2=mpid_2,
            functional=functional,
            external_pressure=external_pressure,
            db_file=db_file,
            high_level=high_level,
        ),
        name="Calculate Adhesion",
        parents=[calc_pes_points, retrieve_matched_slabs],
    )
    fw_list.append(compute_adhesion)

    charge_analysis = Firework(
        FT_StartChargeAnalysisSWF(
            mp_id_1=mpid_1,
            mp_id_2=mpid_2,
            functional=functional,
            external_pressure=external_pressure,
            db_file=db_file,
            high_level=high_level,
        ),
        name="Charge Analysis",
        parents=[calc_pes_points],
    )
    fw_list.append(charge_analysis)

    wf_name = (
        "TriboFlow_"
        + "-".join(
            sorted([f"{formula_1}_({mpid_1})", f"{formula_2}_({mpid_2})"])
        )
        + "_"
        + f"{functional}@{external_pressure}GPa"
    )

    fw_list = Workflow(fw_list, name=wf_name)

    return fw_list
