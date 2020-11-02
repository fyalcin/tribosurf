#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python


from fireworks import Firework
#from triboflow.firetasks.PES import FT_FindHighSymmPoints, FT_StartPESCalcs
from triboflow.firetasks.check_inputs import FT_CheckCompParamDict, \
    FT_CheckInterfaceParamDict, FT_CheckMaterialInputDict, FT_MakeBulkInDB, \
    FT_MakeSlabInDB, FT_MakeInterfaceInDB

__author__ = 'Michael Wolloch'
__copyright__ = 'Copyright 2020, Michael Wolloch'
__version__ = '0.1'
__maintainer__ = 'Michael Wolloch'
__email__ = 'michael.wolloch@univie.ac.at'
__date__ = 'March 11th, 2020'

# =============================================================================
# Custom FireWorks
# =============================================================================

def RunPESCalcsFW(interface_name, functional, tag, FW_name):
    
    FT_1 = FT_FindHighSymmPoints(interface_name=interface_name,
                                 functional=functional)
    
    FT_2 = FT_StartPESCalcs(interface_name=interface_name,
                            functional=functional,
                            tag=tag)
    
    FW = Firework([FT_1, FT_2], name=FW_name)
    
    return FW


def CheckInputsFW(mat1_params, mat2_params, compparams,
                  interface_params, FW_name):
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
    FT_Mat1 = FT_CheckMaterialInputDict(input_dict = mat1_params,
                                        output_dict_name = 'mat_1')
    FT_Mat2 = FT_CheckMaterialInputDict(input_dict = mat2_params,
                                        output_dict_name = 'mat_2')
    
    FT_CompParams = FT_CheckCompParamDict(input_dict = compparams,
                                           output_dict_name = 'comp')
    
    FT_InterParams = FT_CheckInterfaceParamDict(input_dict = interface_params,
                                                output_dict_name = 'inter')
    
    FT_Mat1ToDB = FT_MakeBulkInDB(mat_data_loc = 'mat_1',
                                  comp_data_loc = 'comp')
    FT_Mat2ToDB = FT_MakeBulkInDB(mat_data_loc = 'mat_2',
                                  comp_data_loc = 'comp')
    
    FT_Slab1ToDB = FT_MakeSlabInDB(mat_data_loc = 'mat_1',
                                   comp_data_loc = 'comp')
    FT_Slab2ToDB = FT_MakeSlabInDB(mat_data_loc = 'mat_2',
                                   comp_data_loc = 'comp')
    
    FT_InterfaceToDB = FT_MakeInterfaceInDB(mat1_data_loc = 'mat_1',
                                            mat2_data_loc = 'mat_2',
                                            comp_data_loc = 'comp',
                                            interface_data_loc = 'inter')
    
    FW=Firework([FT_Mat1, FT_Mat2, FT_CompParams, FT_InterParams, FT_Mat1ToDB,
                 FT_Mat2ToDB, FT_Slab1ToDB, FT_Slab2ToDB, FT_InterfaceToDB],
                name = FW_name)
    return FW