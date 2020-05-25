#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python

from fireworks import Firework, ScriptTask
from CommonFiretasks import *
__author__ = 'Michael Wolloch'
__copyright__ = 'Copyright 2020, Michael Wolloch'
__version__ = '0.1'
__maintainer__ = 'Michael Wolloch'
__email__ = 'michael.wolloch@univie.ac.at'
__date__ = 'March 11th, 2020'

# =============================================================================
# Custom FireWorks
# =============================================================================

def ParsePrevVaspCalc_FW(CalcName, OutLoc, PassOnList, PushEnergyLoc=None,
                         spec=None, Name='Parse Vasp Output'):
    """Parse the OUTCAR of a finished VASP caluculation and put it in the spec.
    
    Uses pymatgen to parse the OUTCAR of a previous VASP run and put a
    dictionary containing the final energy, the Fermi energy, the local
    charges, and the local magnetic moments in the fw_spec. If PushEnergyLoc
    is specified, the final energy is also pushed to the end of the array
    at this location in the spec.

    Parameters
    ----------
    CalcName : str
        Label of the previous VASP Calc as given in the CalcLoc.
    OutLoc : list of str
        Each str represents a key specifying the location of the output
        dictionary in the fw_spec.
    PassOnList : list of str
        Specifies the first level keys of the spec that should be passed on to
        the next Firework.
    PushEnergyLoc : list of str, optional
        The parsed total energy will be appended to an array at this location
        in the fw_spec if this is given. The default is None.
    spec : dict, optional
        Previous fw_spec that will be partially updated and/or passed on.
        The default is None.
    Name : str, optional
        Name of the Firework. The default is 'Parse Vasp Output'.

    Returns
    -------
    FW : fireworks.core.firework.Firework
        A Firework parsing VASP OUTCAR output and putting it in the fw_spec.

    """
    from atomate.common.firetasks.glue_tasks import CopyFilesFromCalcLoc
    task_list = []
    task_list.append(CopyFilesFromCalcLoc(calc_loc=CalcName))
    task_list.append(FT_ParseVaspOutput(out_loc=OutLoc,
                                        push_energy_loc=PushEnergyLoc))
    task_list.append(FT_PassSpec(key_list=PassOnList))
    #task_list.append(FT_PrintSpec())
    FW = Firework(task_list, spec = spec, name = Name)
    return FW

def StartDetourWF_FW(WF_type, name, **kwargs):
    """Starts a subworkflow as a detour.
    
    This should be the only way a sub-workflow is started if it should run as
    a detour to the main workflow (inherit children). If a new subworkflow
    is needed, it should be added here as an allowed_type, and the needed
    parameters should be added to the needed_params dictionary.
    

    Parameters
    ----------
    WF_type : str
        Which kind of subworkflow should be executed. Needs to be in
        allowed_types.
    name : str
        Name of the Firework that starts the subworkflow. Not the Subworkflows
        name.
    **kwargs : 
        Pass all the parameters necessary for the subworkflow. Edit
        neede_params if a new subworkflow is added.

    Raises
    ------
    SystemExit
        Will exit if the subworkflow type is not known or the needed_params
        are not passed in **kwargs

    Returns
    -------
    FW : fireworks.core.firework.Firework
        Firework that starts a detour subworkflow via a FWAction.

    """
    # Add more subworkflows here and also the necessary parameters that need
    # to be passed.
    allowed_types = ['Relax_SWF', 'Converge_Kpoints_SWF', 'None']
    needed_params = {'Relax_SWF': ['structure_loc',
                                   'comp_parameters_loc',
                                   'relax_type',
                                   'out_loc',
                                   'to_pass'],
                     'Converge_Kpoints_SWF': ['structure_loc',
                                              'comp_parameters_loc',
                                              'out_loc',
                                              'to_pass'],
                     'None': ['structure_loc', 'out_loc', 'to_pass']}
    if WF_type in allowed_types:
        if not all (k in kwargs for k in needed_params[WF_type]):
            raise SystemExit('Not all necessary arguments passed for the {}'
                             ' Subworkflow.'.format(WF_type))
    else:
        raise SystemExit('Specified subworkflow type {} is not known.\n'
                         'Please select one of: {}'.format(WF_type,
                                                           allowed_types))
    
    if WF_type == 'Relax_SWF':
        FW = Firework(FT_StartRelaxSubWorkflow(
                            structure_loc=kwargs['structure_loc'],
                            comp_parameters_loc=kwargs['comp_parameters_loc'],
                            relax_type=kwargs['relax_type'],
                            out_loc=kwargs['out_loc'],
                            to_pass=kwargs['to_pass']), name=name)
    elif WF_type == 'Converge_Kpoints_SWF':
        FW = Firework(FT_StartKptsConvSubWorkflow(
                            structure_loc=kwargs['structure_loc'],
                            comp_parameters_loc=kwargs['comp_parameters_loc'],
                            out_loc=kwargs['out_loc'],
                            to_pass=kwargs['to_pass']), name=name)
    elif WF_type == 'None':
        FT_copy = FT_CopyInSpec(in_loc=kwargs['structure_loc'],
                                out_loc=kwargs['out_loc'])
        FT_PassOn = FT_PassSpec(key_list = kwargs['to_pass'])
        FW = Firework([FT_copy, FT_PassOn], name=name)
    
    return FW

def TestFW(in_loc, out_loc, pass_list):
    FT0 = FT_PrintSpec()
    FT1 = ScriptTask.from_str('echo "task1"')
    FT2 = FT_CopyInSpec(in_loc=in_loc, out_loc=out_loc)
    FT3 = ScriptTask.from_str('echo "task3"')
    FT4 = ScriptTask.from_str('echo "task4"')
    FT5 = FT_PassSpec(key_list = pass_list)
    FW = Firework([FT1, FT2, FT3, FT5])
    return FW

def FixParametersFW(loc_1, loc_2, out_loc, to_pass, name):
# =============================================================================
#     TODO: Include Encut convergence once it is ready.
# =============================================================================
    FT_K = FT_ChooseCompParams(loc_1=loc_1, loc_2=loc_2, out_loc=out_loc)
    FT_Pass = FT_PassSpec(key_list = to_pass)
    FW = Firework([FT_K, FT_Pass], name=name)
    return FW


def ConvergeParametersFW(key_list, name):
# =============================================================================
#     TODO: Write Firetasks and Fireworks (or better workflows) to converge
#           Kpoints and Energy cutoff.
# =============================================================================
    FW = Firework(FT_PassSpec(key_list=key_list), name=name)
    return FW

def CheckInputsFW(mat1loc, mat2loc, compparamloc, interparamloc, name):
    """Check the input parameters for completeness and add default values.
    
    This Firework uses several Firetasks to check if the necessary input for
    the heterogeneous Triboflow workflow is given and assignes default values
    to optional parameters that are not given. Also makes sure that the
    input parameters are of the correct type and in the correct locations
    in the spec.
    
    Parameters
    ----------
    mat1loc : list of str
        List of keys in the fw_spec pointing to the input dictionary for the
        first material.
    mat2loc : list of str
        List of keys in the fw_spec pointing to the input dictionary for the
        second material.
    compparamloc : list of str
        List of keys in the fw_spec pointing to the input dictionary for the
        computational parameters.
    interparamloc : list of str
        List of keys in the fw_spec pointing to the input dictionary for the
        parameters needed for interface matching. THESE MIGHT CHANGE IF
        INITIALLY NO INTERFACE IS FOUND!
    name : str
        Name for the firework.

    Returns
    -------
    FW : fireworks.core.firework.Firework
        Firework for checking the parameters.

    """
    FT_Mat1 = FT_CheckMaterialInputDict(input_dict_loc=mat1loc)
    FT_Mat2 = FT_CheckMaterialInputDict(input_dict_loc=mat2loc)
    FT_CompParams1 = FT_CheckCompParamDict(input_dict_loc=compparamloc,
                                           output_dict_loc=mat1loc)
    FT_CompParams2 = FT_CheckCompParamDict(input_dict_loc=compparamloc,
                                           output_dict_loc=mat2loc)
    FT_InterParams = FT_CheckInterfaceParamDict(input_dict_loc=interparamloc)
    FT_CompParams3 = FT_CheckCompParamDict(input_dict_loc=compparamloc,
                                output_dict_loc=interparamloc)
    FW=Firework([FT_Mat1, FT_Mat2, FT_CompParams1, FT_CompParams2,
                 FT_InterParams, FT_CompParams3], name=name)
    return FW

    
def ReadWFInputsFW(inputfile_name):
    """Generates a FireWork that reads and checks an input file for a Workflow.
    
    An input file is read in using a FireTask. A following FireTask is used
    to check if essential keys are correctly set in the input file and to
    assign default values for optional keys if they are not given.

    Parameters
    ----------
    inputfile_name : str
        DESCRIPTION.
    essential_keys : list of str, optional
        DESCRIPTION.
    optional_keys : list of str, optional
        DESCRIPTION.

    Returns
    -------
    FireWork
    """
    
    import os
    indir = os.getcwd()
    fw = Firework([FT_ReadInputFile(filename=inputfile_name),
               FT_CheckInputDict(input_dict_name=inputfile_name)],
              {'_launch_dir': indir})
    
    return fw


def MakeGeneralizedKpointsFW(structure, PrecalcDict, IncarDict=None,
                              WorkDir='./', KPTS_loc=None):
    """
    A workflow is constructured that uses the K-Point Grid Generator
    of the Mueller group at John Hopkins to get a kpoint mesh for
    a given structure. Details can be found at
    http://muellergroup.jhu.edu/K-Points.html

    Args:
        structure (pymatgen structure object): Required input structure for
                                               which the Kpoints mesh should
                                               be created. 
        PrecalcDict (dict): Required parameters for the Grid Generator, should
                            definately include the MINDISTANCE flag.
        IncarDict (dict):   Optional INCAR file as a dictionary to pass
                            information about ISYM, MAGMOM, etc...
        WorkDir (str):      Optional working directory where everything happens
                            this is created if not already found.
        KPTS_loc (list of str): Optional location of the Kpoints object in
                                the spec. Will default to 'KPOINTS' if left at
                                'None'.

    Returns:
        Firework
    """
    import os

    if WorkDir != './':
        if not os.path.isdir(WorkDir):
            os.mkdir(WorkDir)
        
    
    FT_write_files = FT_WriteGeneralKpointsInputs(structure=structure,
                                                  PrecalcDict=PrecalcDict,
                                                  IncarDict=IncarDict,
                                                  WorkDir=WorkDir)
    FT_makeKpoints = ScriptTask(script=['(cd '+WorkDir+' ; getKPoints)'])
    FT_getKpoints = FT_GetKpoints(WorkDir=WorkDir, KPTS_loc=KPTS_loc)
    FT_CleanUp = FT_CleanUpGeneralizedKpointFiles(WorkDir=WorkDir)

    FW = Firework([FT_write_files, FT_makeKpoints, FT_getKpoints,
                   FT_CleanUp])

    return FW



def MakeHeteroInterfaceFromScratchFW(Material_1, Material_2, parameters):
    """
    A FireWork is constructed that forms two slabs out of information
    contained in Material_1 and Material_2. These slabs are then combined
    to a heterogeneous interface. Functions from MPInterfaces are used.
    If no match can be found, the parameters handling the lattice matching
    in the parameters dictionary are gradually increased by 5% until a match
    can be made. An output structure is not only saved to the database but
    also moved to the input directory.

    Args:
        Material_1 (dict):  Dictionary containing at least the keys (values):
                            formula (str), miller (list), thickness (float).
                            Might also contain MP_ID (str) and vacuum (float).
                            E.g.:
                                    {'formula': 'Fe2O3',
                                     'miller': [1, 1, 0],
                                     'thickness': 12.5,
                                     'MP_ID': 'mp-19770',
                                     'vacuum': 20.0}
        Material_2 (dict):  Dictionary containing at least the keys (values):
                            formula (str), miller (list), thickness (float).
                            Might also contain MP_ID (str) and vacuum (float).
                            E.g.:
                                    {'formula': 'C',
                                     'miller': [1, 1, 1],
                                     'thickness': 12.5,
                                     'MP_ID': 'mp-66',
                                     'vacuum': 20.0}
        parameters (dict):  Dictionary containing the key (values):
                            interface_distance (float): Distance between the
                                two materials.
                            max_area (float): Maximal cross section
                                area of the matched cell in Angstrom squared.
                                Defaults to 200.,
                            max_mismatch (float): Maximal allowed mismatch
                                between lattice vector length in %.
                                Defaults to 0.01,
                            max_angle_diff (float): Maximal allowed mismatch
                                between lattice vector angle in degrees.
                                Defaults to 1.
                            r1r2_tol (float): Tolerance between factors r1 and
                                r2 which relate to the matched cell areas as:
                                abs(r1*area1-r2*area2)/max_area <= r1r2_tol.
                                Defaults to 0.02
    Returns:
        Firework
    """
    #get the current directory to later copy the output structure to here.
    
    #Make the bottom slab
    formula = Material_1['formula']
    miller = Material_1['miller']
    thickness = Material_1['thickness']
    if 'MP_ID' in Material_1:
        mp_id = Material_1['MP_ID']
    else:
        mp_id = None
    if 'vacuum' in Material_1:
        vacuum = Material_1['vacuum']
    else:
        vacuum = None
    
    miller_compact = ''.join(str(e) for e in miller)
    bottom_slab_name = formula + miller_compact
    FT_bottom_slab = FT_MakeSlabFromFormula(formula=formula, miller=miller,
                                            thickness=thickness, MP_ID=mp_id,
                                            vacuum=vacuum)
    
    #Make the top slab
    formula = Material_2['formula']
    miller = Material_2['miller']
    thickness = Material_2['thickness']
    if 'MP_ID' in Material_2:
        mp_id = Material_2['MP_ID']
    else:
        mp_id = None
    if 'vacuum' in Material_2:
        vacuum = Material_2['vacuum']
    else:
        vacuum = None
        
    miller_compact = ''.join(str(e) for e in miller)
    top_slab_name = formula + miller_compact
    FT_top_slab = FT_MakeSlabFromFormula(formula=formula, miller=miller,
                                         thickness=thickness, MP_ID=mp_id,
                                         vacuum=vacuum)    
    
    #Make the matched interface
    params={'interface_distance': parameters['interface_distance']}
    
    if 'max_area' in parameters:
        params['max_area'] = parameters['max_area']
    else:
        params['max_area'] = 200.0
        
    if 'max_mismatch' in parameters:
        params['max_mismatch'] = parameters['max_mismatch']
    else:
        params['max_mismatch'] = 0.01
    
    if 'max_angle_diff' in parameters:
        params['max_angle_diff'] = parameters['max_angle_diff']
    else:
        params['max_angle_diff'] = 1.0
        
    if 'r1r2_tol' in parameters:
        params['r1r2_tol'] = parameters['r1r2_tol']
    else:
        params['r1r2_tol'] = 0.02
        
    FT_get_hetero_structure = FT_MakeHeteroStructure(
                                bottom_slab_name=bottom_slab_name,
                                top_slab_name=top_slab_name,
                                parameters=params)
    
    import os
    indir = os.getcwd()
    
    FW = Firework([FT_bottom_slab, FT_top_slab, FT_get_hetero_structure],
                  {'_launch_dir': indir})
    
    return FW