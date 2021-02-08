#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python


from fireworks import Firework, FileTransferTask, ScriptTask
from triboflow.firetasks.PES import FT_FindHighSymmPoints, FT_StartPESCalcs, \
    FT_RetrievePESEnergies, FT_ComputePES
from triboflow.firetasks.check_inputs import FT_CheckCompParamDict, \
    FT_CheckInterfaceParamDict, FT_CheckMaterialInputDict, FT_MakeBulkInDB, \
    FT_MakeSlabInDB, FT_MakeInterfaceInDB
from triboflow.utils.file_manipulation import copy_output_files

__author__ = 'Michael Wolloch'
__copyright__ = 'Copyright 2020, Michael Wolloch'
__version__ = '0.1'
__maintainer__ = 'Michael Wolloch'
__email__ = 'michael.wolloch@univie.ac.at'
__date__ = 'March 11th, 2020'

# =============================================================================
# Custom FireWorks
# =============================================================================

def run_pes_calc_fw(top_slab, bottom_slab, interface_name,
                    functional, comp_parameters, tag, FW_name):
    """Compute high-symmetry points for an interface and start PES calculations.
    
    Combines two Fireworks that find the high-symmetry points for the interface
    and start the VASP calculations for the unique high-symmetry points
    respectively.
    

    Parameters
    ----------
    top_slab : pymatgen.core.surface.Slab
        Top slab of the interface.
    bottom_slab : pymatgen.core.surface.Slab
        Bottom slab of the interface.
    interface_name : str
        Unique name for the interface that is used in the output and the
        database.
    functional : str
        Functional to be used. 'PBE' or 'SCAN' will work.
    comp_parameters : dict
        Dictionary containing computational options. E.g. encut, k_dens, vdw,...
    tag : str
        combination of the interface_name and a uuid. To uniquely identify the
        computations in the database.
    FW_name : str
        Name of the Firework.

    Returns
    -------
    FW : fireworks.core.firework.Firework
        First Firework of a PES subworkflow.

    """
    FT_1 = FT_FindHighSymmPoints(top_slab=top_slab, bot_slab=bottom_slab,
                                 interface_name=interface_name,
                                 functional=functional)
    
    FT_2 = FT_StartPESCalcs(top_slab=top_slab, bot_slab=bottom_slab,
                            interface_name=interface_name,
                            functional=functional,
                            comp_parameters=comp_parameters,
                            tag=tag)
    
    FW = Firework([FT_1, FT_2], name=FW_name)
    
    return FW

def make_pes_fw(interface_name, functional, tag, FW_name, file_output,
                output_dir, remote_copy=False, server=None, user=None, 
                port=None):
    """Retrieve PES calculations from the database and compute the PES.
    
    Retrive the computed energies of the unique high-symmetry points and match
    them to the replicate points. Duplicates the points, interpolates with
    radial basis functions and saves the results. Plots the results as well.
    Optionally write file output and copy it to a output directory.
    

    Parameters
    ----------
    interface_name : str
        Unique name for the interface that is used in the output and the
        database.
    functional : str
        Functional to be used. 'PBE' or 'SCAN' will work.
    tag : str
        combination of the interface_name and a uuid. To uniquely identify the
        computations in the database.
    FW_name : str
        Name of the Firework.
    file_output : bool
        Determines if files are written to disk.
    output_dir : str
        Location the output files are copied to if file_output is selected.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The user name on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.
        
    Returns
    -------
    FW : fireworks.core.firework.Firework
        Final Firework of a PES subworkflow.

    """    
    FT_1 = FT_RetrievePESEnergies(interface_name=interface_name,
                                  functional=functional,
                                  tag=tag)
    FT_2 = FT_ComputePES(interface_name=interface_name,
                         functional=functional,
                         file_output=file_output)
    
    if file_output: 
        output_files = ['Computet_PES_data_'+interface_name+'.dat',
                        'Interpolated_PES_data_'+interface_name+'.dat',
                        'PES_' + str(interface_name) + '.png']
        FT_3 = copy_output_files(file_list=output_files,
                                 output_dir=output_dir,
                                 remote_copy=remote_copy,
                                 server=server,
                                 user=server,
                                 port=port)
            
        FW = Firework([FT_1, FT_2, FT_3], name=FW_name)
    else:
        
        FW = Firework([FT_1, FT_2], name=FW_name)
    
    return FW


def check_inputs_fw(mat1_params, mat2_params, compparams,
                    interface_params, FW_name):
    """Check the input parameters for completeness and add default values.
    
    This Firework uses several Firetasks to check if the necessary input for
    the heterogeneous Triboflow workflow is given and assigns default values
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
    FT_Mat1 = FT_CheckMaterialInputDict(input_dict=mat1_params,
                                        output_dict_name='mat_1')
    FT_Mat2 = FT_CheckMaterialInputDict(input_dict=mat2_params,
                                        output_dict_name='mat_2')
    
    FT_CompParams = FT_CheckCompParamDict(input_dict=compparams,
                                          output_dict_name='comp')
    
    FT_InterParams = FT_CheckInterfaceParamDict(input_dict=interface_params,
                                                output_dict_name='inter')
    
    FT_Mat1ToDB = FT_MakeBulkInDB(mat_data_loc='mat_1',
                                  comp_data_loc='comp')
    FT_Mat2ToDB = FT_MakeBulkInDB(mat_data_loc='mat_2',
                                  comp_data_loc='comp')
    
    FT_Slab1ToDB = FT_MakeSlabInDB(mat_data_loc='mat_1',
                                   comp_data_loc='comp')
    FT_Slab2ToDB = FT_MakeSlabInDB(mat_data_loc='mat_2',
                                   comp_data_loc='comp')
    
    FT_InterfaceToDB = FT_MakeInterfaceInDB(mat1_data_loc='mat_1',
                                            mat2_data_loc='mat_2',
                                            comp_data_loc='comp',
                                            interface_data_loc='inter')
    
    FW=Firework([FT_Mat1, FT_Mat2, FT_CompParams, FT_InterParams, FT_Mat1ToDB,
                 FT_Mat2ToDB, FT_Slab1ToDB, FT_Slab2ToDB, FT_InterfaceToDB],
                name=FW_name)
    
    return FW