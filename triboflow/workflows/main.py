"""Main Workflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from fireworks import Workflow, Firework
from triboflow.common_firetasks import FT_AddSelectiveDynamics, FT_WritePrecalc, \
    FT_PrintSpec, FT_StructFromVaspOutput, FT_FetchStructureFromFormula, \
    FT_MakeSlabFromStructure, FT_MakeHeteroStructure, FT_PassSpec, \
    FT_LoopKpoints, FT_PassKpointsInfo, FT_EnergyCutoffConvo

from triboflow.helper_functions import GetCustomVaspRelaxSettings, \
    GetCustomVaspStaticSettings

def Heterogeneous_WF(inputs):
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
    WF = []

    Initialize = Firework(FT_PrintSpec(), spec=inputs,
                                    name= 'Initialize Workflow')
    WF.append(Initialize)
    
    Check_Inputs = CheckInputsFW(mat1loc=['material_1'], mat2loc=['material_2'],
                                 compparamloc=['computational_params'],
                                 interparamloc=['interface_params'],
                                 name='Check Inputs')
    WF.append(Check_Inputs)
    
    Start_M1 = Firework(FT_FetchStructureFromFormula(
                            materials_dict_loc=['material_1']),
                            name='Get '+inputs['material_1']['formula']+' bulk')
    WF.append(Start_M1)
    
    Start_M2 = Firework(FT_FetchStructureFromFormula(
                            materials_dict_loc=['material_2']),
                            name='Get '+inputs['material_2']['formula']+' bulk')
    WF.append(Start_M2)
    
    bulk_loc = ['material_1', inputs['material_1']['formula']+'_fromMP']
    ConvergeEncut_M1 = StartDetourWF_FW('Converge_Encut_SWF',
                            name='Converge Encut '+
                            inputs['material_1']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_1'],
                            out_loc=['material_1'],
                            to_pass=['material_1'])
    WF.append(ConvergeEncut_M1)
    
    bulk_loc = ['material_1', inputs['material_1']['formula']+'_equiVol']
    Converge_M1 = StartDetourWF_FW('Converge_Kpoints_SWF',
                            name='Converge Kpoints '+
                            inputs['material_1']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_1'],
                            out_loc=['material_1'],
                            to_pass=['material_1'])
    WF.append(Converge_M1)
    
    bulk_loc = ['material_2', inputs['material_2']['formula']+'_fromMP']
    ConvergeEncut_M2 = StartDetourWF_FW('Converge_Encut_SWF',
                            name='Converge Encut '+
                            inputs['material_2']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_2'],
                            out_loc=['material_2'],
                            to_pass=['material_2', 'computational_params',
                                     'interface_params'])
    WF.append(ConvergeEncut_M2)
    
    bulk_loc = ['material_2', inputs['material_2']['formula']+'_equiVol']
    Converge_M2 = StartDetourWF_FW('Converge_Kpoints_SWF',
                            name='Converge Kpoints '+
                            inputs['material_2']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_2'],
                            out_loc=['material_2'],
                            to_pass=['material_2', 'computational_params',
                                     'interface_params'])
    WF.append(Converge_M2)
    
    Final_Params = FixParametersFW(loc_1=['material_1'],
                                   loc_2=['material_2'],
                                   out_loc=['interface_params'],
                                   to_pass=['_all'],
                                   name='Select computational parameters')
    WF.append(Final_Params)
    
    bulk_loc = ['material_1', inputs['material_1']['formula']+'_equiVol']
    out_loc = ['material_1', inputs['material_1']['formula']+'_relaxed']
    Relax_M1 = StartDetourWF_FW('None',
                            name='Relax '+inputs['material_1']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_1'],
                            relax_type='bulk_full_relax',
                            out_loc=out_loc,
                            to_pass=['material_1', 'computational_params',
                                     'interface_params'])
    WF.append(Relax_M1)
    
    bulk_loc = ['material_2', inputs['material_2']['formula']+'_equiVol']
    out_loc = ['material_2', inputs['material_2']['formula']+'_relaxed']
    Relax_M2 = StartDetourWF_FW('None',
                            name='Relax '+inputs['material_2']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_2'],
                            relax_type='bulk_full_relax',
                            out_loc=out_loc,
                            to_pass=['material_2', 'computational_params',
                                     'interface_params'])
    WF.append(Relax_M2)
    
    bulk_loc = ['material_1', inputs['material_1']['formula']+'_relaxed']
    Make_Slab_M1 = Firework(FT_MakeSlabFromStructure(bulk_loc=bulk_loc,
                                                     dict_loc=['material_1']),
                    name='Make '+inputs['material_1']['formula']+' slab')
    WF.append(Make_Slab_M1)
    
    bulk_loc = ['material_2', inputs['material_2']['formula']+'_relaxed']
    Make_Slab_M2 = Firework(FT_MakeSlabFromStructure(bulk_loc=bulk_loc,
                                                     dict_loc=['material_2']),
                    name='Make '+inputs['material_2']['formula']+' slab')
    WF.append(Make_Slab_M2)
    
    bottom_slab_loc = ['material_1', inputs['material_1']['formula']+
                         inputs['material_1']['miller']]
    out_loc = ['material_1', inputs['material_1']['formula']+
                         inputs['material_1']['miller']+'_relaxed']
    Relax_Slab_M1 = StartDetourWF_FW('Relax_SWF',
                            name='Relax '+inputs['material_1']['formula']+
                            ' slab',
                            structure_loc=bottom_slab_loc,
                            comp_parameters_loc=['material_1'],
                            relax_type='slab_pos_relax',
                            out_loc=out_loc,
                            to_pass=['material_1', 'computational_params',
                                     'interface_params'])
    WF.append(Relax_Slab_M1)
    
    top_slab_loc = ['material_2', inputs['material_2']['formula']+
                         inputs['material_2']['miller']]
    out_loc = ['material_2', inputs['material_2']['formula']+
                         inputs['material_2']['miller']+'_relaxed']
    Relax_Slab_M2 = StartDetourWF_FW('Relax_SWF',
                            name='Relax '+inputs['material_2']['formula']+
                            ' slab',
                            structure_loc=top_slab_loc,
                            comp_parameters_loc=['material_2'],
                            relax_type='slab_pos_relax',
                            out_loc=out_loc,
                            to_pass=['material_2', 'computational_params',
                                     'interface_params'])
    WF.append(Relax_Slab_M2)
    
    bottom_slab_loc = ['material_1', inputs['material_1']['formula']+
                         inputs['material_1']['miller']+'_relaxed']
    top_slab_loc = ['material_2', inputs['material_2']['formula']+
                         inputs['material_2']['miller']+'_relaxed']
    Make_Hetero_Structure = Firework(FT_MakeHeteroStructure(
        bottom_slab_loc=bottom_slab_loc, top_slab_loc=top_slab_loc,
        parameters_loc=['interface_params']), name='Make the interface')
    WF.append(Make_Hetero_Structure)
    
    Print_spec = Firework(FT_PrintSpec(), name='Print Spec at the end')
    WF.append(Print_spec)
    
    #Define dependencies:
    Dependencies = {Initialize: [Check_Inputs],
                    Check_Inputs: [Start_M1, Start_M2],
                    Start_M1: [ConvergeEncut_M1],
                    Start_M2: [ConvergeEncut_M2],
                    ConvergeEncut_M1: [Converge_M1],
                    ConvergeEncut_M2: [Converge_M2],
                    Converge_M1: [Final_Params],
                    Converge_M2: [Final_Params],
                    Final_Params: [Relax_M1, Relax_M2],
                    Relax_M1: [Make_Slab_M1],
                    Relax_M2: [Make_Slab_M2],
                    Make_Slab_M1: [Relax_Slab_M1],
                    Make_Slab_M2: [Relax_Slab_M2],
                    Relax_Slab_M1: [Make_Hetero_Structure],
                    Relax_Slab_M2: [Make_Hetero_Structure],
                    Make_Hetero_Structure: [Print_spec]}

    WF = Workflow(WF, Dependencies, name='Heterogeneous Workflow')
    return WF

