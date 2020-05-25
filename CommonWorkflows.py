#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Common Workflows to be used by the FireFlow Project.
Created on Thu Apr 16 15:46:22 2020

@author: mwo
"""


from fireworks import Workflow, ScriptTask, Firework, FileTransferTask
from atomate.vasp.config import VASP_CMD
from atomate.utils.utils import env_chk
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.common.firetasks.glue_tasks import PassCalcLocs, \
    CopyFilesFromCalcLoc
from CommonFiretasks import FT_AddSelectiveDynamics, FT_WritePrecalc, \
    FT_PrintSpec, FT_StructFromVaspOutput, FT_FetchStructureFromFormula, \
    FT_MakeSlabFromStructure, FT_MakeHeteroStructure, FT_PassSpec, \
    FT_LoopKpoints, FT_PassKpointsInfo
from CommonFireworks import CheckInputsFW, ConvergeParametersFW, \
    FixParametersFW, StartDetourWF_FW, ParsePrevVaspCalc_FW
from HelperFunctions import GetCustomVaspRelaxSettings, \
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
    Converge_M1 = StartDetourWF_FW('Converge_Kpoints_SWF',
                            name='Converge Kpoints '+
                            inputs['material_1']['formula'],
                            structure_loc=bulk_loc,
                            comp_parameters_loc=['material_1'],
                            out_loc=['material_1'],
                            to_pass=['material_1'])
    WF.append(Converge_M1)
    
    bulk_loc = ['material_2', inputs['material_2']['formula']+'_fromMP']
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
    
    bulk_loc = ['material_1', inputs['material_1']['formula']+'_fromMP']
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
    
    bulk_loc = ['material_2', inputs['material_2']['formula']+'_fromMP']
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
    Relax_Slab_M1 = StartDetourWF_FW('None',
                            name='Relax '+inputs['material_1']['formula']+
                            ' slab',
                            structure_loc=bottom_slab_loc,
                            comp_parameters_loc=['material_2'],
                            relax_type='slab_pos_relax',
                            out_loc=out_loc,
                            to_pass=['material_1', 'computational_params',
                                     'interface_params'])
    WF.append(Relax_Slab_M1)
    
    top_slab_loc = ['material_2', inputs['material_2']['formula']+
                         inputs['material_2']['miller']]
    out_loc = ['material_2', inputs['material_2']['formula']+
                         inputs['material_2']['miller']+'_relaxed']
    Relax_Slab_M2 = StartDetourWF_FW('None',
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
                    Start_M1: [Converge_M1],
                    Start_M2: [Converge_M2],
                    Converge_M2: [Final_Params],
                    Converge_M1: [Final_Params],
                    Final_Params: [Relax_M1, Relax_M2],
                    Relax_M1: [Make_Slab_M1],
                    Relax_M2: [Make_Slab_M2],
                    Make_Slab_M1: [Relax_Slab_M1],
                    Make_Slab_M2: [Relax_Slab_M2],
                    Relax_Slab_M1: [Make_Hetero_Structure],
                    Relax_Slab_M2: [Make_Hetero_Structure],
                    Make_Hetero_Structure: [Print_spec]}

    WF = Workflow(WF, Dependencies, name='Dummy Heterogeneous Workflow')
    return WF

def ConvergeKpoints_SWF(structure, comp_parameters, out_loc, to_pass, spec,
                        k_dist_incr=1.0, n_converge=3):
    """Subworkflows that converges the k_distance for generalized MP-meshes.
    
    Takes a given structure and computational parameters and makes consecutive
    VASP calculations with increasingly dense generalized Monkhorst-Pack
    K-meshes until the total energy is converged. Uses an external script and
    K-point grid generator from the Mueller group at John Hopkins
    http://muellergroup.jhu.edu/K-Points.html

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-pint grids.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
    out_loc : list of str
        DESCRIPTION.
    to_pass : list of str
        List of keys that each represent a location in the first level of the
        fw_spec and signifies which of those are to be passed to the next
        Firework. E.g. if the fw_sec = {'a': a, 'b': {'b1'; b1, 'b2': b2},
        'c': c}} and to_pass = ['a', 'b'], the spec in the next FW will be:
        {'a': a, 'b': {'b1'; b1, 'b2': b2}}
    spec : dict
        Previous fw_spec that will be updated and/or passed on for child
        Fireworks.
    k_dist_incr : float, optional
        Increment for the k_distance during the convergence. Defaults to 1.0
    n_converge : int, optional
        Number of calculations that have to show the same energy as the last
        one as to signify convergence, Defaults to 3.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A subworkflow intended to find the converged k_distance for a given
        structure.

    """
    
    name = 'Kpoint_Convergence SWF of '+structure.composition.reduced_formula
    
    if 'energy_tolerance' not in comp_parameters:
        comp_parameters['energy_tolerance'] = 0.001
    
    FT_Loop_Kpoints = FT_LoopKpoints(structure = structure,
                                   out_loc = out_loc,
                                   comp_params = comp_parameters,
                                   k_dist_incr = k_dist_incr,
                                   n_converge = n_converge)
    FT_PassOn = FT_PassSpec(key_list=to_pass)
    
    FW = Firework([FT_Loop_Kpoints, FT_PassOn], spec=spec,
                  name = 'Start Convergence')
    FW_2 = Firework([FT_PassKpointsInfo(out_loc=out_loc), FT_PassOn],
                    spec=spec, name = 'Prepare Output and Pass Spec')
    WF = Workflow([FW, FW_2], {FW: [FW_2]}, name=name)
    return WF
    
def GetEnergy_SWF(structure, comp_parameters, static_type, out_loc,
                  to_pass, spec, push_energy_loc=None):
    """Return a subworkflow that computes the total energy of a given structure.
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-pint grids.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
    static_type : str
        Specifies in which way the vasp input parameters are set. Check
        HelperFunction.GetCustomVaspStaticSettings for details.
    out_loc : list of str
        List of keys that specify in which location the output data will be
        stored.
    to_pass : list of str
        List of keys that each represent a location in the first level of the
        fw_spec and signifies which of those are to be passed to the next
        Firework. E.g. if the fw_sec = {'a': a, 'b': {'b1'; b1, 'b2': b2},
        'c': c}} and to_pass = ['a', 'b'], the spec in the next FW will be:
        {'a': a, 'b': {'b1'; b1, 'b2': b2}}
    spec : dict
        Previous fw_spec that will be updated and/or passed on for child
        Fireworks.
    push_energy_loc : list of str, optional
        The parsed total energy will be appended to an array at this location
        in the fw_spec if this is given. The default is None.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A subworkflow intended to compute the total energy for a given
        structure using VASP.

    """
    CalcName = structure.composition.reduced_formula+' GetEnergy'
    
    vis, uis, vdw = GetCustomVaspStaticSettings(structure, comp_parameters,
                                                static_type)
            
    custom_params = {'user_incar_settings': uis,
                     'vdw': vdw,
                     'user_potcar_functional': 'PBE_54'}
    
    FT_list = []
    
    FT_list.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vis,
                             vasp_input_params=custom_params))
    
    if 'k_distance' in comp_parameters:
        Precalc_Dict = {'MINDISTANCE': comp_parameters['k_distance'],
                        'MINTOTALKPOINTS': 4,
                        'GAPDISTANCE': 6,
                        'MONOCLINIC_SEARCH_DEPTH': 2500,
                        'TRICLINIC_SEARCH_DEPTH': 1500}
        FT_list.append(FT_WritePrecalc(PrecalcDict=Precalc_Dict))
        FT_list.append(ScriptTask(script=['getKPoints']))
    elif comp_parameters.get('functional') == 'SCAN':
        if 'is_metal' in comp_parameters:
            if comp_parameters.get('is_metal'):
                k_distance = 50
            else:
                k_distance = 31.5
        else:
            k_distance = 50
        Precalc_Dict = {'MINDISTANCE': k_distance,
                        'MINTOTALKPOINTS': 4,
                        'GAPDISTANCE': 6,
                        'MONOCLINIC_SEARCH_DEPTH': 2500,
                        'TRICLINIC_SEARCH_DEPTH': 1500}
        FT_list.append(FT_WritePrecalc(PrecalcDict=Precalc_Dict))
        FT_list.append(ScriptTask(script=['getKPoints']))
        
    if vdw:
        vdw_kernel = env_chk('>>vdw_kernel_dir<<', spec)+'/vdw_kernel.bindat'
        FT_list.append(FileTransferTask({'files': [{'src': vdw_kernel,
                                            'dest': 'vdw_kernel.bindat'}],
                                         'mode': 'copy'}))
    
    FT_list.append(RunVaspCustodian(vasp_cmd=VASP_CMD, auto_npar=True))
    
    FT_list.append(VaspToDb(db_file='>>db_file<<'))
    FT_list.append(PassCalcLocs(name=CalcName))
    
    FW_Static = Firework(FT_list, name=CalcName+' FW')
    
    FW_PP = ParsePrevVaspCalc_FW(CalcName=CalcName, OutLoc=out_loc,
                                 PassOnList=to_pass,
                                 PushEnergyLoc=push_energy_loc,
                                 Name='Parse Output FW', spec=spec)
    
    WF = Workflow([FW_Static, FW_PP], {FW_Static: [FW_PP]}, name=CalcName+
                                                              ' Workflow')
    return WF

def Relax_SWF(structure, comp_parameters, relax_type, out_loc, to_pass, spec):
    """Relax bulk, slab, or interface structures in a subworkflow.
    
    A FireWorks (Sub)workflow is created that works with the PBE and SCAN
    functionals. It also supports generalized Kpoint-grids and different
    kinds of settings depending on the 'relax_type' input parameters. For
    Interfaces, selective dynamics will be added so that ions can only move
    in the direciton of the third lattice. Functionality can be added to pass
    a custom list for selective dynamics as well. The workflow puts the final
    relaxed structure in the spec of the firework which is then passed on to
    the main workflow.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure to relax.
    comp_parameters : dict
        Dictionary of computational parameters for the relaxation.
    relax_type : str
        Used by the 'GetCustomVaspRelaxSettings' HelperFunction to set
        override parameters for the used MPRelaxSer or the MPScanRelaxSet and
        also determine which form of relaxation (ISIF) to do.
    out_loc : list of str
        List of keys that specify in which location the final output structure
        will be stored in the spec.
    to_pass : list of str
        List of all the keys in the spec to be passed to the next firework.
        May be ['_all'] to pass the whole spec.
    spec : dict
        present spec that will be passed on to further fireworks after the
        relaxed structure has been added.

    Returns
    -------
    WF : FireWorks Workflow

    """  
    CalcName = structure.composition.reduced_formula+relax_type
    
    vis, uis, vdw = GetCustomVaspRelaxSettings(structure, comp_parameters,
                                               relax_type)
    
            
    custom_params = {'user_incar_settings': uis,
                     'vdw': vdw,
                     'user_potcar_functional': 'PBE_54'}
    
    FT_list = []
    
    FT_list.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vis,
                             vasp_input_params=custom_params))
    
    if relax_type.startswith('interface'):
        FT_list = FT_AddSelectiveDynamics(structure=structure)
    
    if 'k_distance' in comp_parameters:
        Precalc_Dict = {'MINDISTANCE': comp_parameters['k_distance'],
                        'MINTOTALKPOINTS': 4,
                        'GAPDISTANCE': 6,
                        'MONOCLINIC_SEARCH_DEPTH': 2500,
                        'TRICLINIC_SEARCH_DEPTH': 1500}
        FT_list.append(FT_WritePrecalc(PrecalcDict=Precalc_Dict))
        FT_list.append(ScriptTask(script=['getKPoints']))
    elif comp_parameters.get('functional') == 'SCAN':
        if 'is_metal' in comp_parameters:
            if comp_parameters.get('is_metal'):
                k_distance = 50
            else:
                k_distance = 31.5
        else:
            k_distance = 50
        Precalc_Dict = {'MINDISTANCE': k_distance,
                        'MINTOTALKPOINTS': 4,
                        'GAPDISTANCE': 6,
                        'MONOCLINIC_SEARCH_DEPTH': 2500,
                        'TRICLINIC_SEARCH_DEPTH': 1500}
        FT_list.append(FT_WritePrecalc(PrecalcDict=Precalc_Dict))
        FT_list.append(ScriptTask(script=['getKPoints']))
    
    if 'functional' in comp_parameters:
        if comp_parameters['functional'] == 'SCAN':
            job_type = 'metagga_opt_run'
        else:
            job_type = 'double_relaxation_run'
    else:
        job_type = 'double_relaxation_run'
        
    if vdw:
        vdw_kernel = env_chk('>>vdw_kernel_dir<<', spec)+'/vdw_kernel.bindat'
        FT_list.append(FileTransferTask({'files': [{'src': vdw_kernel,
                                            'dest': 'vdw_kernel.bindat'}],
                                         'mode': 'copy'}))
    
    FT_list.append(RunVaspCustodian(vasp_cmd=VASP_CMD, job_type=job_type,
                                  ediffg=uis['EDIFFG'],
                                  auto_npar=True,
                                  half_kpts_first_relax=False))
    
    FT_list.append(VaspToDb(db_file='>>db_file<<'))
    FT_list.append(PassCalcLocs(name=CalcName))
    
    FW_Relax = Firework(FT_list, name=CalcName+' FW')
    
    
    FW_PP = Firework([CopyFilesFromCalcLoc(calc_loc=CalcName,
                                           filenames=['OUTCAR.relax2.gz',
                                                      'CONTCAR.relax2.gz']),
                      FT_StructFromVaspOutput(out_struct_loc=out_loc),
                      FT_PassSpec(key_list=to_pass)],
                     spec=spec, name=CalcName+' post processing FW')
    
    WF = Workflow([FW_Relax, FW_PP], {FW_Relax: [FW_PP]}, name=CalcName+
                                                              ' Workflow')
    return WF
    
    
     
     
    