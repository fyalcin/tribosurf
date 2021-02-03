#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 15:41:22 2021

Modules to put nulk, slabs, and interfaces in DB during initializing phase.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__credits__ = 'This module is based on the Triboflow package, Michael Wolloch'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'January 20th, 2021'

from triboflow.utils.database import Navigator, NavigatorMP
from triboflow.tasks.handle_struct import interface_name

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase
from atomate.utils.utils import env_chk


# ============================================================================
# FireTasks
# ============================================================================

@explicit_serialize
class FTPutMaterialInDB(FiretaskBase):
    """ 
    Put a bulk and the required slab from MP into the high level DB.
    
    Parameters
    ----------
    mat : str
        Location of the dictionary with the material data in the spec.
    comp_paramas : str
        Location of the dictionary with computational parameters in the spec.
    db_file : str
        Location of the high level db on the working machine.
    
    """
    
    _fw_name = 'Make bulk and interface entry into high level DB'
    required_params = ['mat', 'comp_params']
    optional_params = ['db_file']
    
    def run_task(self, fw_spec):
        """ Run the FireTask.
        """

        # Extract data from input dictionaries
        data = fw_spec[self['mat']]
        comp_params = fw_spec[self['comp_params']]
        
        # Locate the high level DB
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        # Put bulk and slab in high level DB
        put_bulk_in_db(data, comp_params, db_file)
        put_slab_in_db(data, comp_params, db_file)

@explicit_serialize
class FTPutInterfaceInDB(FiretaskBase):
    """
    Put the information to make an interface between two materials into the 
    high level DB.
    
    Parameters
    ----------
    mat_1 : str
        Location of the dictionary with the first material data in the spec.
    mat_2 : str
        Location of the dictionary with the second material data in the spec.
    comp_paramas : str
        Location of the dictionary with computational parameters in the spec.
    inter_paramas : str
        Location of the dictionary with interfacial parameters in the spec.
    db_file : str
        Location of the high level db on the working machine.
    
    """
    
    _fw_name = 'Make interface entry into high level DB'
    required_params = ['mat_1', 'mat_2', 'comp_params', 'inter_params']
    optional_params = ['db_file']
    
    def run_task(self, fw_spec):
        """ Run the FireTask.
        """

        # Extract data from input dictionaries
        data_1 = fw_spec[self['mat_1']]
        data_2 = fw_spec[self['mat_2']]
        comp_params = fw_spec[self['comp_params']]
        inter_params = fw_spec[self['inter_params']]

        # Locate the high level DB        
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        # Put interface info and parameters in high level DB
        put_inter_in_db(data_1, data_2, comp_params, inter_params, db_file)
    
    
# ============================================================================
# Function to put the bulk and slab of a given materials in the user DB
# ============================================================================

def put_bulk_in_db(data, comp_params, db_file):
    """
    Put bulk data in high level DB.

    Parameters
    ----------
    data : dict
        Dictionary with bulk data.
    comp_paramas : str
        Computational parameters of the calculation.
    db_file : str
        Local DataBase on the machine.

    """
    
    # Extract data from dictionaries and from MP using NavigatorMP
    struct, mp_id, functional, comp_params = material_data_dor_db(
        data, comp_params, db_file, metal_thr = 0.2)
    
    # Load the bulk structure into the high level DB
    nav = Navigator(db_file, high_level=True)
    nav.insert_data(functional+'.bulk_data', 
                    {'mpid': mp_id,
                     'formula': data['formula'],
                     'structure_fromMP': struct.as_dict(),
                     'comp_parameters': comp_params})
    
def put_slab_in_db(data, comp_params, db_file):
    """
    Put slab data in high level DB.

    Parameters
    ----------
    data : dict
        Dictionary with slab data.
    comp_paramas : str
        Computational parameters of the calculation.
    db_file : str
        Local DataBase on the machine.

    """

    # Extract data from dictionaries and from MP using NavigatorMP
    _, mp_id, functional, comp_params = material_data_dor_db(
        data, comp_params, db_file, metal_thr = 0.5)

    # Load the slab structure into the high level DB
    nav = Navigator(db_file, high_level=True)
    nav.insert_data(functional+'.slab_data', 
                    {'mpid': mp_id,
                     'formula': data['formula'],
                     'miller': data['miller'],
                     'min_thickness': data['min_thickness'],
                     'min_vacuum': data['min_vacuum'],
                     'comp_parameters': comp_params})

def put_inter_in_db(data_1, data_2, comp_params, inter_params, db_file):
    """
    Put interfacial data in high level DB.

    Parameters
    ----------
    data_1: dict
        Dictionary with the first material data.
    data_2: dict
        Dictionary with the second material data.
    comp_paramas: str
        Computational parameters of the calculation.
    inter_params: str
        Interfacial parameters.
    db_file: str
        Local DataBase on the machine.

    """

    # Extract data from dictionaries and from MP using NavigatorMP
    functional = comp_params['functional']
    struct_1, mp_id_1, _, _ = material_data_dor_db(data_1, comp_params, 
                                                   db_file, metal_thr=None)
    struct_2, mp_id_2, _, _ = material_data_dor_db(data_2, comp_params,
                                                   db_file, metal_thr=None)
    name = interface_name(mp_id_1, data_1['miller'], mp_id_2, data_2['miller'])
    
    # Load the interface structure into the high level DB
    nav_high = Navigator(db_file, high_level=True)
    nav_high.insert_data(functional+'.interface_data', 
                         {'name': name,
                          'comp_parameters': comp_params,
                          'interface_parameters': inter_params})

def material_data_dor_db(data, comp_params, db_file, metal_thr=None):
    """
    Extract material information from data dictionary and MP online database.
    Returned parameters are suitable to be stored in high level DB.

    Parameters
    ----------
    data : dict
        Dictionary containing the data about the material.
    comp_paramas : str
        Computational parameters of the calculation.
    db_file : str
        Local DataBase on the machine.
    metal_thr : float or None, optional
        Float threshold to be compared to the bandgap of the material, to 
        understand if it is a metal or not. The default is None.

    Returns
    -------
    struct : pymatgen.Structure
        Material structure downloaded from Materials Project DB.
    mp_id : str
        MP ID for the selected material.
    functional : str
        Functional of the calculation.
    comp_params : dict
        Computational parameters for the calculation, it is different from the
        input value only if metal_thr is set.

    """
    
    # Extract relevant data
    mp_id = data['mpid']
    formula = data['formula']
    functional = comp_params['functional']
    
    # Check an online bandgap and understand if material is metal or not
    nav_mp = NavigatorMP()
    struct, mp_id = nav_mp.get_low_energy_structure(formula, mp_id)
    
    # TODO : Is it correct that you need to slice the returned value?
    bandgap = nav_mp.get_property_from_mp(mp_id, properties=['band_gap'])['band_gap']
    
    # Set the metal value for the given material
    if metal_thr is not None:
        comp_params['is_metal'] = True if bandgap <= metal_thr else False
    
    return struct, mp_id, functional, comp_params
