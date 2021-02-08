"""Main Workflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from fireworks import Workflow, Firework

from triboflow.fireworks.common import check_inputs_fw
from triboflow.firetasks.encut_convergence import FT_StartEncutConvo
from triboflow.firetasks.kpoint_convergence import FT_StartKPointConvo
from triboflow.firetasks.structure_manipulation import FT_StartSlabRelaxSWF, \
    FT_MakeHeteroStructure
from triboflow.firetasks.PES import FT_StartPESCalcSubWF
from triboflow.firetasks.check_inputs import FT_UpdateCompParams
from triboflow.utils.database import NavigatorMP
from triboflow.utils.structure_manipulation import interface_name

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
    mat_1 = inputs.get('material_1')
    mat_2 = inputs.get('material_2')
    comp_params = inputs.get('computational_params')
    inter_params = inputs.get('interface_params')
    
    if not all([mat_1, mat_2, comp_params, inter_params]):
        raise SystemExit('The inputs-dictionary for this workflow must '
                         'contain the following keys:\n'
                         'material_1\nmaterial_2\ncomputational_params\n'
                         'interface_params')

    nav_mp = NavigatorMP()
    struct_1, mp_id_1 = nav_mp.get_low_energy_structure(
        chem_formula=mat_1.get('formula'),
        mp_id=mat_1.get('mp_id'))
    struct_2, mp_id_2 = nav_mp.get_low_energy_structure(
        chem_formula=mat_2.get('formula'),
        mp_id=mat_2.get('mp_id'))
    
    functional = comp_params.get('functional', 'SCAN')
    
    WF = []

    Initialize = check_inputs_fw(mat1_params=mat_1,
                                 mat2_params=mat_2,
                                 compparams=comp_params,
                                 interface_params=inter_params,
                                 FW_name='Check input parameters')
    WF.append(Initialize)
    
    ConvergeEncut_M1 = Firework(FT_StartEncutConvo(mp_id=mp_id_1,
                                                   functional=functional),
                                name='Start encut convergence for {}'
                                      .format(mat_1['formula']))
    WF.append(ConvergeEncut_M1)
    
    ConvergeEncut_M2 = Firework(FT_StartEncutConvo(mp_id=mp_id_2,
                                                   functional=functional),
                                name='Start encut convergence for {}'
                                      .format(mat_2['formula']))
    WF.append(ConvergeEncut_M2)
    
    ConvergeKpoints_M1 = Firework(FT_StartKPointConvo(mp_id=mp_id_1,
                                                   functional=functional),
                                  name='Start kpoints convergence for {}'
                                        .format(mat_1['formula']))
    WF.append(ConvergeKpoints_M1)
    
    ConvergeKpoints_M2 = Firework(FT_StartKPointConvo(mp_id=mp_id_2,
                                                      functional=functional),
                                name='Start kpoints convergence for {}'
                                     .format(mat_2['formula']))
    WF.append(ConvergeKpoints_M2)
    
    Final_Params = Firework(FT_UpdateCompParams(mp_id_1=mp_id_1,
                                                mp_id_2=mp_id_2,
                                                miller_1=mat_1.get('miller'),
                                                miller_2=mat_2.get('miller'),
                                                functional=functional),
                            name='Consolidate computational paramters')
    WF.append(Final_Params)
    
    MakeSlabs_M1 = Firework(FT_StartSlabRelaxSWF(mp_id=mp_id_1,
                                                 miller=mat_1.get('miller'),
                                                 functional=functional),
                            name='Make and relax slabs for {}'
                                 .format(mat_1['formula']))
    WF.append(MakeSlabs_M1)
    
    MakeSlabs_M2 = Firework(FT_StartSlabRelaxSWF(mp_id=mp_id_2,
                                                 miller=mat_2.get('miller'),
                                                 functional=functional),
                            name='Make and relax slabs for {}'
                                  .format(mat_2['formula']))
    WF.append(MakeSlabs_M2)
    
    MakeInterface = Firework(FT_MakeHeteroStructure(
        mp_id_1=mp_id_1,
        mp_id_2=mp_id_2,
        miller_1=mat_1.get('miller'),
        miller_2=mat_2.get('miller'),
        functional=functional),
        name='Match the interface')
    WF.append(MakeInterface)
    
    CalcPESPoints = Firework(FT_StartPESCalcSubWF(
        mp_id_1=mp_id_1,
        mp_id_2=mp_id_2,
        miller_1=mat_1.get('miller'),
        miller_2=mat_2.get('miller'),
        functional=functional),
        name='Compute PES high-symmetry points')
    WF.append(CalcPESPoints)
    
    # Define dependencies:
    Dependencies = {Initialize: [ConvergeEncut_M1, ConvergeEncut_M2],
                    ConvergeEncut_M1: [ConvergeKpoints_M1],
                    ConvergeEncut_M2: [ConvergeKpoints_M2],
                    ConvergeKpoints_M1: [Final_Params],
                    ConvergeKpoints_M2: [Final_Params],
                    Final_Params: [MakeSlabs_M1, MakeSlabs_M2],
                    MakeSlabs_M1: [MakeInterface],
                    MakeSlabs_M2: [MakeInterface],
                    MakeInterface: [CalcPESPoints]}

    WF_Name = 'TriboFlow ' + interface_name(mp_id_1, mat_1.get('miller'),
                                            mp_id_2, mat_2.get('miller'))

    WF = Workflow(WF, Dependencies, name=WF_Name)

    return WF