"""Main Workflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from fireworks import Workflow, Firework

from triboflow.fireworks.init_fws import InitWF
from triboflow.firetasks.structure_manipulation import (
    FT_MakeHeteroStructure,
    FT_StartBulkPreRelax,
)
from triboflow.firetasks.init_check import unbundle_input, material_from_mp
from triboflow.firetasks.utils import FT_CopyHomogeneousSlabs
from triboflow.firetasks.check_inputs import FT_UpdateInterfaceCompParams
from triboflow.firetasks.adhesion import (
    FT_RelaxMatchedSlabs,
    FT_RetrieveMatchedSlabs,
)
from triboflow.firetasks.start_swfs import (
    FT_StartAdhesionSWF,
    FT_StartBulkConvoSWF,
    FT_StartDielectricSWF,
    FT_StartPESCalcSWF,
    FT_StartChargeAnalysisSWF,
)
from triboflow.utils.structure_manipulation import interface_name
from triboflow.firetasks.run_slabs_wfs import FT_SlabOptThick


def homogeneous_wf(inputs):
    mat, comp_params, interface_params = unbundle_input(
        inputs, keys=["material", "computational_params", "interface_params"]
    )
    struct, mp_id = material_from_mp(mat)
    functional = comp_params.get("functional", "PBE")

    WF = []

    # define max area here large, even if it will not matter since this is a
    # homogeneous interface. However, max_area is a necessary key.

    # pressure might default to None, so we have to check for that
    pressure = interface_params.get("external_pressure", 0.0) or 0.0

    Initialize = InitWF.checkinp_homo_interface(
        material=mat,
        computational=comp_params,
        interface={"max_area": 200, "external_pressure": pressure},
    )
    WF.append(Initialize)

    PreRelaxation = Firework(
        FT_StartBulkPreRelax(mp_id=mp_id, functional=functional),
        name=f'Start pre-relaxation for {mat["formula"]}',
    )
    WF.append(PreRelaxation)

    ConvergeEncut = Firework(
        FT_StartBulkConvoSWF(
            conv_type="encut", mp_id=mp_id, functional=functional
        ),
        name=f'Start encut convergence for {mat["formula"]}',
    )
    WF.append(ConvergeEncut)

    ConvergeKpoints = Firework(
        FT_StartBulkConvoSWF(
            conv_type="kpoints", mp_id=mp_id, functional=functional
        ),
        name=f'Start kpoints convergence for {mat["formula"]}',
    )
    WF.append(ConvergeKpoints)

    CalcDielectric = Firework(
        FT_StartDielectricSWF(
            mp_id=mp_id,
            functional=functional,
            update_bulk=True,
            update_slabs=True,
        ),
        name=f'Start dielectric SWF for {mat["formula"]}',
    )
    #WF.append(CalcDielectric)

    Final_Params = Firework(
        FT_UpdateInterfaceCompParams(
            mp_id_1=mp_id,
            mp_id_2=mp_id,
            miller_1=mat.get("miller"),
            miller_2=mat.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Consolidate computational parameters",
    )
    WF.append(Final_Params)

    MakeSlabs = Firework(
        FT_SlabOptThick(
            mp_id=mp_id, miller=mat.get("miller"), functional=functional
        ),
        name=f'Slab thickness optimization for {mat["formula"]}',
    )
    WF.append(MakeSlabs)

    MakeInterface = Firework(
        FT_MakeHeteroStructure(
            mp_id_1=mp_id,
            mp_id_2=mp_id,
            miller_1=mat.get("miller"),
            miller_2=mat.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Match the interface",
    )
    WF.append(MakeInterface)

    CopySlabs = Firework(
        FT_CopyHomogeneousSlabs(
            mp_id=mp_id,
            miller=mat.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Copy homogeneous slabs",
    )
    WF.append(CopySlabs)

    CalcPESPoints = Firework(
        FT_StartPESCalcSWF(
            mp_id_1=mp_id,
            mp_id_2=mp_id,
            miller_1=mat.get("miller"),
            miller_2=mat.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Compute PES high-symmetry points",
    )
    WF.append(CalcPESPoints)

    ComputeAdhesion = Firework(
        FT_StartAdhesionSWF(
            mp_id_1=mp_id,
            mp_id_2=mp_id,
            miller_1=mat.get("miller"),
            miller_2=mat.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Calculate Adhesion",
    )
    WF.append(ComputeAdhesion)

    ChargeAnalysis = Firework(
        FT_StartChargeAnalysisSWF(
            mp_id_1=mp_id,
            mp_id_2=mp_id,
            miller_1=mat.get("miller"),
            miller_2=mat.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Charge Analysis",
    )
    WF.append(ChargeAnalysis)

    # Define dependencies:
    Dependencies = {
        Initialize: [PreRelaxation],
        PreRelaxation: [ConvergeEncut],
        ConvergeEncut: [ConvergeKpoints],
        #ConvergeKpoints: [CalcDielectric],
        #CalcDielectric: [Final_Params],
        ConvergeKpoints: [Final_Params],
        Final_Params: [MakeSlabs],
        MakeSlabs: [MakeInterface],
        MakeInterface: [CalcPESPoints, CopySlabs],
        CalcPESPoints: [ComputeAdhesion, ChargeAnalysis],
        CopySlabs: [ComputeAdhesion],
    }

    WF_Name = (
        "TriboFlow_"
        + interface_name(mp_id, mat.get("miller"), mp_id, mat.get("miller"))
        + "_"
        + f"{functional}@{pressure}GPa"
    )

    WF = Workflow(WF, Dependencies, name=WF_Name)

    return WF


def heterogeneous_wf(inputs):
    """Return main workflow for heterogeneous interfaces within Triboflow.

    Parameters
    ----------
    inputs : dict
        Dictionary containing sub-dictionaries with material information and
        parameters for the interface creation and computational settings.

    Returns
    -------
    WF : FireWorks Workflow
        Main Triboflow workflow for heterogeneous interfaces.

    """
    mat_1, mat_2, comp_params, inter_params = unbundle_input(inputs)

    struct_1, mp_id_1 = material_from_mp(mat_1)
    struct_2, mp_id_2 = material_from_mp(mat_2)

    functional = comp_params.get("functional", "PBE")

    WF = []

    # pressure might default to None, so we have to check for that
    pressure = inter_params.get("external_pressure", 0.0) or 0.0

    Initialize = InitWF.checkinp_hetero_interface(
        material_1=mat_1,
        material_2=mat_2,
        computational=comp_params,
        interface=inter_params,
    )
    WF.append(Initialize)

    PreRelaxation_M1 = Firework(
        FT_StartBulkPreRelax(mp_id=mp_id_1, functional=functional),
        name="Start pre-relaxation for {}".format(mat_1["formula"]),
    )
    WF.append(PreRelaxation_M1)

    PreRelaxation_M2 = Firework(
        FT_StartBulkPreRelax(mp_id=mp_id_2, functional=functional),
        name="Start pre-relaxation for {}".format(mat_2["formula"]),
    )
    WF.append(PreRelaxation_M2)

    ConvergeEncut_M1 = Firework(
        FT_StartBulkConvoSWF(
            conv_type="encut",
            mp_id=mp_id_1,
            functional=functional,
        ),
        name="Start encut convergence for {}".format(mat_1["formula"]),
    )
    WF.append(ConvergeEncut_M1)

    ConvergeEncut_M2 = Firework(
        FT_StartBulkConvoSWF(
            conv_type="encut",
            mp_id=mp_id_2,
            functional=functional,
        ),
        name="Start encut convergence for {}".format(mat_2["formula"]),
    )
    WF.append(ConvergeEncut_M2)

    ConvergeKpoints_M1 = Firework(
        FT_StartBulkConvoSWF(
            conv_type="kpoints",
            mp_id=mp_id_1,
            functional=functional,
        ),
        name="Start kpoints convergence for {}".format(mat_1["formula"]),
    )
    WF.append(ConvergeKpoints_M1)

    ConvergeKpoints_M2 = Firework(
        FT_StartBulkConvoSWF(
            conv_type="kpoints",
            mp_id=mp_id_2,
            functional=functional,
        ),
        name="Start kpoints convergence for {}".format(mat_2["formula"]),
    )
    WF.append(ConvergeKpoints_M2)

    CalcDielectric_M1 = Firework(
        FT_StartDielectricSWF(
            mp_id=mp_id_1,
            functional=functional,
            update_bulk=True,
            update_slabs=True,
        ),
        name=f'Start dielectric SWF for {mat_1["formula"]}',
    )
    #WF.append(CalcDielectric_M1)

    CalcDielectric_M2 = Firework(
        FT_StartDielectricSWF(
            mp_id=mp_id_2,
            functional=functional,
            update_bulk=True,
            update_slabs=True,
        ),
        name=f'Start dielectric SWF for {mat_2["formula"]}',
    )
    #WF.append(CalcDielectric_M2)

    Final_Params = Firework(
        FT_UpdateInterfaceCompParams(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get("miller"),
            miller_2=mat_2.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Consolidate computational parameters",
    )
    WF.append(Final_Params)

    MakeSlabs_M1 = Firework(
        FT_SlabOptThick(
            mp_id=mp_id_1, miller=mat_1.get("miller"), functional=functional
        ),
        name="Slab thickness optimization " "for {}".format(mat_1["formula"]),
    )
    WF.append(MakeSlabs_M1)

    MakeSlabs_M2 = Firework(
        FT_SlabOptThick(
            mp_id=mp_id_2, miller=mat_2.get("miller"), functional=functional
        ),
        name="Slab thickness optimization " "for {}".format(mat_2["formula"]),
    )
    WF.append(MakeSlabs_M2)

    MakeInterface = Firework(
        FT_MakeHeteroStructure(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get("miller"),
            miller_2=mat_2.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Match the interface",
    )
    WF.append(MakeInterface)

    RelaxMatchedSlabs = Firework(
        FT_RelaxMatchedSlabs(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get("miller"),
            miller_2=mat_2.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Fully relax the matched slabs",
    )
    WF.append(RelaxMatchedSlabs)

    RetrieveMatchedSlabs = Firework(
        FT_RetrieveMatchedSlabs(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get("miller"),
            miller_2=mat_2.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Retrieve relaxed matched slabs",
    )
    WF.append(RetrieveMatchedSlabs)

    CalcPESPoints = Firework(
        FT_StartPESCalcSWF(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get("miller"),
            miller_2=mat_2.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Compute PES high-symmetry points",
    )
    WF.append(CalcPESPoints)

    ComputeAdhesion = Firework(
        FT_StartAdhesionSWF(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get("miller"),
            miller_2=mat_2.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Calculate Adhesion",
    )
    WF.append(ComputeAdhesion)

    ChargeAnalysis = Firework(
        FT_StartChargeAnalysisSWF(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get("miller"),
            miller_2=mat_2.get("miller"),
            functional=functional,
            external_pressure=pressure,
        ),
        name="Charge Analysis",
    )
    WF.append(ChargeAnalysis)

    # Define dependencies:
    Dependencies = {
        Initialize: [PreRelaxation_M1, PreRelaxation_M2],
        PreRelaxation_M1: [ConvergeEncut_M1],
        PreRelaxation_M2: [ConvergeEncut_M2],
        ConvergeEncut_M1: [ConvergeKpoints_M1],
        ConvergeEncut_M2: [ConvergeKpoints_M2],
        # ConvergeKpoints_M1: [CalcDielectric_M1],
        # ConvergeKpoints_M2: [CalcDielectric_M2],
        # CalcDielectric_M1: [Final_Params],
        # CalcDielectric_M2: [Final_Params],
        ConvergeKpoints_M1: [Final_Params],
        ConvergeKpoints_M2: [Final_Params],
        Final_Params: [MakeSlabs_M1, MakeSlabs_M2],
        MakeSlabs_M1: [MakeInterface],
        MakeSlabs_M2: [MakeInterface],
        MakeInterface: [CalcPESPoints, RelaxMatchedSlabs],
        RelaxMatchedSlabs: [RetrieveMatchedSlabs],
        CalcPESPoints: [ComputeAdhesion, ChargeAnalysis],
        RetrieveMatchedSlabs: [ComputeAdhesion],
    }

    WF_Name = (
        "TriboFlow_"
        + interface_name(
            mp_id_1, mat_1.get("miller"), mp_id_2, mat_2.get("miller")
        )
        + "_"
        + f"{functional}@{pressure}GPa"
    )

    WF = Workflow(WF, Dependencies, name=WF_Name)

    return WF
