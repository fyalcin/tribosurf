#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Common Workflows to be used by the FireFlow Project.
Created on Thu Apr 16 15:46:22 2020

@author: mwo
"""


from fireworks import Workflow, ScriptTask, Firework, FileTransferTask
from atomate.vasp.config import RELAX_MAX_FORCE, VASP_CMD, DB_FILE
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.common.firetasks.glue_tasks import PassCalcLocs, \
    CopyFilesFromCalcLoc
#from custodian.vasp.jobs import double_relaxation_run, metagga_opt_run
from CommonFireworks import FT_AddSelectiveDynamics, FT_WritePrecalc, \
    FT_PrintSpec, FT_StructFromVaspOutput
from HelperFunctions import GetCustomVaspRelaxSettings


def Relax_SWF(structure, comp_parameters, relax_type, out_loc, spec):
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
    spec : dict
        present spec that will be passed on to further fireworks after the
        relaxed structure has been added.

    Returns
    -------
    WF : FireWorks Workflow

    """  
    CalcName = structure.composition.reduced_formula+relax_type
    
    vis, uis, vdw = GetCustomVaspRelaxSettings(comp_parameters, relax_type)
            
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
    
    if 'functional' in comp_parameters:
        if comp_parameters['functional'] == 'SCAN':
            job_type = 'metagga_opt_run'
        else:
            job_type = 'double_relaxation_run'
    else:
        job_type = 'double_relaxation_run'
        
    if vdw:
        vdw_kernel = '/home/mwo/bin/vdw_kernel/vdw_kernel.bindat'
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
                      FT_PrintSpec()],
                     spec=spec, name=CalcName+' post processing FW')
    
    WF = Workflow([FW_Relax, FW_PP], {FW_Relax: [FW_PP]}, name=CalcName+
                                                              ' Workflow')
    return WF
    
    
     
     
    