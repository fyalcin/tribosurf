"""Main Workflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from fireworks import Workflow, Firework
from triboflow.fireworks.common import CheckInputsFW
from triboflow.firetasks.encut_convergence import FT_StartEncutConvo
from triboflow.firetasks.kpoint_convergence import FT_StartKPointConvo
from triboflow.helper_functions import GetLowEnergyStructure

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
    mat_1 = inputs.get('material_1')
    mat_2 = inputs.get('material_2')
    comp_params = inputs.get('computational_params')
    inter_params = inputs.get('interface_params')
    
    if not all([mat_1, mat_2, comp_params, inter_params]):
        raise SystemExit('The inputs-dictionary for this workflow must '
                         'contain the following keys:\n'
                         'material_1\nmaterial_2\ncomputational_params\n'
                         'interface_params')
        
    struct_1, mp_id_1 = GetLowEnergyStructure(mat_1.get('formula'),
                                              mat_1.get('mp_id'))
    struct_2, mp_id_2 = GetLowEnergyStructure(mat_1.get('formula'),
                                              mat_1.get('mp_id'))
    
    functional = comp_params.get('functional', 'SCAN')
    
    WF = []

    Initialize = CheckInputsFW(mat1_params = mat_1,
                               mat2_params = mat_2,
                               compparams = comp_params,
                               interface_params = inter_params,
                               FW_name = 'Check input parameters FW')
    WF.append(Initialize)
    
    ConvergeEncut_M1 = Firework(FT_StartEncutConvo(mp_id = mp_id_1,
                                                   functional = functional),
                                name = 'Start encut convergence for {}'
                                        .format(mat_1['formula']))
    WF.append(ConvergeEncut_M1)
    
    ConvergeEncut_M2 = Firework(FT_StartEncutConvo(mp_id = mp_id_2,
                                                   functional = functional),
                                name = 'Start encut convergence for {}'
                                        .format(mat_2['formula']))
    WF.append(ConvergeEncut_M2)
    
    ConvergeKpoints_M1 = Firework(FT_StartKPointConvo(mp_id = mp_id_1,
                                                   functional = functional),
                                name = 'Start kpoints convergence for {}'
                                        .format(mat_1['formula']))
    WF.append(ConvergeKpoints_M1)
    
    ConvergeKpoints_M2 = Firework(FT_StartKPointConvo(mp_id = mp_id_2,
                                                   functional = functional),
                                name = 'Start kpoints convergence for {}'
                                        .format(mat_2['formula']))
    WF.append(ConvergeKpoints_M2)
    
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

