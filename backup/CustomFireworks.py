#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python

from fireworks import Firework, ScriptTask
from fireworks.core.firework import FWAction, FiretaskBase
from CustomFiretasks import *
__author__ = 'Michael Wolloch'
__copyright__ = 'Copyright 2020, Michael Wolloch'
__version__ = '0.1'
__maintainer__ = 'Michael Wolloch'
__email__ = 'michael.wolloch@univie.ac.at'
__date__ = 'March 11th, 2020'

# =============================================================================
# Custom FireWorks
# =============================================================================

def MakeGeneralizedKpoints_FW(structure, PrecalcDict, IncarDict=None,
                              WorkDir='./'):
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

    Returns:
        Firework
    """
    import os

    if WorkDir != './':
        if not os.path.isdir(WorkDir):
            os.mkdir(WorkDir)
        
    
    FT_write_files = WriteGeneralKpointsInputs(structure=structure,
                                               PrecalcDict=PrecalcDict,
                                               IncarDict=IncarDict,
                                               WorkDir=WorkDir)
    FT_makeKpoints = ScriptTask(script=['(cd '+WorkDir+' ; getKPoints)'])
    FT_getKpoints = GetKpoints(WorkDir=WorkDir)
    FT_CleanUp = CleanUpGeneralizedKpointFiles(WorkDir=WorkDir)

    FW = Firework([FT_write_files, FT_makeKpoints, FT_getKpoints,
                   FT_CleanUp])

    return FW


def MakeHeteroInterface(Material_1, Material_2, parameters):
    """
    A FireWork is constructed that forms two slabs out of information
    contained in Material_1 and Material_2. These slabs are then combined
    to a heterogeneous interface. Functions from MPInterfaces are used.

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
    FT_bottom_slab = MakeSlabFromFormula(formula=formula, miller=miller,
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
    FT_top_slab = MakeSlabFromFormula(formula=formula, miller=miller,
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
        
    FT_get_hetero_structure = MakeHeteroStructure(
                                bottom_slab_name=bottom_slab_name,
                                top_slab_name=top_slab_name,
                                parameters=params)
    
    
    FW = Firework([FT_bottom_slab, FT_top_slab, FT_get_hetero_structure])
    
    return FW