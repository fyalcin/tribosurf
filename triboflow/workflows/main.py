"""Main Workflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from fireworks import Workflow, Firework

from triboflow.fireworks.init_fws import InitWF
from triboflow.firetasks.structure_manipulation import (
    FT_MakeHeteroStructure, FT_StartBulkPreRelax)
from triboflow.firetasks.init_check import unbundle_input, material_from_mp
from triboflow.firetasks.check_inputs import FT_UpdateInterfaceCompParams
from triboflow.firetasks.adhesion import (
    FT_RelaxMatchedSlabs, FT_RetrieveMatchedSlabs)
from triboflow.firetasks.start_swfs import (FT_StartAdhesionSWF,
                                            FT_StartBulkConvoSWF,
                                            FT_StartDielectricSWF,
                                            FT_StartPESCalcSWF,
                                            FT_StartChargeAnalysisSWF)
from triboflow.utils.structure_manipulation import interface_name
from triboflow.firetasks.run_slabs_wfs import FT_SlabOptThick

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
    
    functional = comp_params.get('functional', 'PBE')
    
    WF = []

    Initialize = InitWF.checkinp_hetero_interface(material_1=mat_1,
                                                  material_2=mat_2,
                                                  computational=comp_params,
                                                  interface=inter_params)
    WF.append(Initialize)

    PreRelaxation_M1 = Firework(FT_StartBulkPreRelax(mp_id=mp_id_1,
                                                     functional=functional),
                                name='Start pre-relaxation for {}'
                                .format(mat_1['formula']))
    WF.append(PreRelaxation_M1)

    PreRelaxation_M2 = Firework(FT_StartBulkPreRelax(mp_id=mp_id_2,
                                                     functional=functional),
                                name='Start pre-relaxation for {}'
                                .format(mat_2['formula']))
    WF.append(PreRelaxation_M2)

    ConvergeEncut_M1 = Firework(FT_StartBulkConvoSWF(conv_type='encut',
                                                     mp_id=mp_id_1,
                                                     functional=functional,
                                                     ),
                                name='Start encut convergence for {}'
                                .format(mat_1['formula']))
    WF.append(ConvergeEncut_M1)
    
    ConvergeEncut_M2 = Firework(FT_StartBulkConvoSWF(conv_type='encut',
                                                     mp_id=mp_id_2,
                                                     functional=functional,
                                                     ),
                                name='Start encut convergence for {}'
                                .format(mat_2['formula']))
    WF.append(ConvergeEncut_M2)
    
    ConvergeKpoints_M1 = Firework(FT_StartBulkConvoSWF(conv_type='kpoints',
                                                       mp_id=mp_id_1,
                                                       functional=functional,
                                                       ),
                                  name='Start kpoints convergence for {}'
                                  .format(mat_1['formula']))
    WF.append(ConvergeKpoints_M1)
    
    ConvergeKpoints_M2 = Firework(FT_StartBulkConvoSWF(conv_type='kpoints',
                                                       mp_id=mp_id_2,
                                                       functional=functional,
                                                       ),
                                  name='Start kpoints convergence for {}'
                                  .format(mat_2['formula']))
    WF.append(ConvergeKpoints_M2)
    
    CalcDielectric_M1 = Firework(FT_StartDielectricSWF(mp_id=mp_id_1,
                                                       functional=functional,
                                                       update_bulk=True,
                                                       update_slabs=True),
                                 name=f'Start dielectric SWF for {mat_1["formula"]}')
    WF.append(CalcDielectric_M1)
    
    CalcDielectric_M2 = Firework(FT_StartDielectricSWF(mp_id=mp_id_2,
                                                       functional=functional,
                                                       update_bulk=True,
                                                       update_slabs=True),
                                 name=f'Start dielectric SWF for {mat_2["formula"]}')
    WF.append(CalcDielectric_M2)
    
    
    Final_Params = Firework(FT_UpdateInterfaceCompParams(mp_id_1=mp_id_1,
                                                         mp_id_2=mp_id_2,
                                                         miller_1=mat_1.get('miller'),
                                                         miller_2=mat_2.get('miller'),
                                                         functional=functional),
                            name='Consolidate computational parameters')
    WF.append(Final_Params)
    
    MakeSlabs_M1 = Firework(FT_SlabOptThick(mp_id=mp_id_1,
                                            miller=mat_1.get('miller'),
                                            functional=functional),
                                            name='Slab thickness optimization '
                                                 'for {}'.format(mat_1['formula']))
    WF.append(MakeSlabs_M1)

    MakeSlabs_M2 = Firework(FT_SlabOptThick(mp_id=mp_id_2,
                                            miller=mat_2.get('miller'),
                                            functional=functional),
                                            name='Slab thickness optimization '
                                                 'for {}'.format(mat_2['formula']))
    WF.append(MakeSlabs_M2)
    
    MakeInterface = Firework(FT_MakeHeteroStructure(
        mp_id_1=mp_id_1,
        mp_id_2=mp_id_2,
        miller_1=mat_1.get('miller'),
        miller_2=mat_2.get('miller'),
        functional=functional),
        name='Match the interface')
    WF.append(MakeInterface)
    
    RelaxMatchedSlabs = Firework(
        FT_RelaxMatchedSlabs(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get('miller'),
            miller_2=mat_2.get('miller'),
            functional=functional),
        name='Fully relax the matched slabs')
    WF.append(RelaxMatchedSlabs)
    
    RetrieveMatchedSlabs = Firework(
        FT_RetrieveMatchedSlabs(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get('miller'),
            miller_2=mat_2.get('miller'),
            functional=functional),
        name='Retrieve relaxed matched slabs')
    WF.append(RetrieveMatchedSlabs)
    
    CalcPESPoints = Firework(
        FT_StartPESCalcSWF(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get('miller'),
            miller_2=mat_2.get('miller'),
            functional=functional),
        name='Compute PES high-symmetry points')
    WF.append(CalcPESPoints)
    
    ComputeAdhesion = Firework(
        FT_StartAdhesionSWF(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get('miller'),
            miller_2=mat_2.get('miller'),
            functional=functional),
        name='Calculate Adhesion')
    WF.append(ComputeAdhesion)
    
    ChargeAnalysis = Firework(
        FT_StartChargeAnalysisSWF(
            mp_id_1=mp_id_1,
            mp_id_2=mp_id_2,
            miller_1=mat_1.get('miller'),
            miller_2=mat_2.get('miller'),
            functional=functional),
        name='Charge Analysis')
    WF.append(ChargeAnalysis)
    
    # Define dependencies:
    Dependencies = {Initialize: [PreRelaxation_M1, PreRelaxation_M2],
                    PreRelaxation_M1: [ConvergeEncut_M1],
                    PreRelaxation_M2: [ConvergeEncut_M2],
                    ConvergeEncut_M1: [ConvergeKpoints_M1],
                    ConvergeEncut_M2: [ConvergeKpoints_M2],
                    ConvergeKpoints_M1: [CalcDielectric_M1],
                    ConvergeKpoints_M2: [CalcDielectric_M2],
                    CalcDielectric_M1: [Final_Params],
                    CalcDielectric_M2: [Final_Params],
                    Final_Params: [MakeSlabs_M1, MakeSlabs_M2],
                    MakeSlabs_M1: [MakeInterface],
                    MakeSlabs_M2: [MakeInterface],
                    MakeInterface: [CalcPESPoints, RelaxMatchedSlabs],
                    RelaxMatchedSlabs: [RetrieveMatchedSlabs],
                    CalcPESPoints: [ComputeAdhesion, ChargeAnalysis],
                    RetrieveMatchedSlabs: [ComputeAdhesion]}

    WF_Name = 'TriboFlow ' + interface_name(mp_id_1, mat_1.get('miller'),
                              mp_id_2, mat_2.get('miller')) +' '+functional

    WF = Workflow(WF, Dependencies, name=WF_Name)

    return WF
