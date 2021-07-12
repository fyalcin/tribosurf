#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:35:43 2021

This modules contains classes and functions useful to initialize different
prebuilt and/or custom workflows.

The module contains:

** InitWF **:
    Collection of static methods tot check the input parameters for a Workflow.
    It includes the following methods:
        - checkinp_hetero_interface
        - checkinp_homo_interface
    
    Author: Gabriele Losi (glosi000)
    Credits: The code is partially based on the original work of Michael 
    Wolloch, Triboflow package, Wien University
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = 'Gabriele Losi'
__credits__ = 'This module is based on the work of Michael Wolloch, TriboFlow'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'January 20th, 2021'

from fireworks import Firework
from triboflow.firetasks.init_check import FTCheckInput
from triboflow.firetasks.init_db import (
    FT_PutMaterialInDB, 
    FT_PutInterfaceInDB
)

class InitWF:
    """
    This is a collection of static methods that check the input parameters for 
    a Workflow and check their completeness, adding default values when needed.
    
    """
    
    @staticmethod
    def checkinp_hetero_interface(material_1, material_2, computational,
                                  interface, fw_name='Check input parameters'):        
        """
        Create a Fireworks to check if the necessary input for a heterogenous
        workflow are given and assignes default values to optional parameters. 
        Input parameters are checked for correct type and location in the spec.
        
        Parameters
        ----------
        material_1 : list of str
            Keys list in fw_spec pointing to the input dict of first material.

        material_2 : list of str
            Keys list in fw_spec pointing to the input dict of second material.

        computational : list of str
            Keys list in fw_spec pointing to the computational parameters.

        interface : list of str
            Keys list in fw_spec pointing to the interface parameters, which
            are needed for interface matching.

        fw_name : str
            Name of the returned FireWork.
    
        Returns
        -------
        FW : fireworks.core.firework.Firework
            Firework that checks all input parameters for an heterogeneous WF.
            
        """

        # Firetasks checking the materials parameters
        ft_mat1 = FTCheckInput(input_dict = material_1,
                               read_key = 'material_params',
                               output_dict_name = 'mat_1')
        ft_mat2 = FTCheckInput(input_dict = material_2,
                               read_key = 'material_params',
                               output_dict_name = 'mat_2')
        
        # Firetask checking the computational parameters
        ft_computation = FTCheckInput(input_dict = computational, 
                                      read_key = 'comp_params',
                                      output_dict_name = 'comp')
        
        # Firetask checking the interfacial matching parameters
        ft_interface = FTCheckInput(input_dict = interface,
                                    read_key = 'interface_params',
                                    output_dict_name = 'inter')
        
        # Put materials bulk and slab in DB
        ft_mat1_db = FT_PutMaterialInDB(mat = 'mat_1', comp_params = 'comp')
        ft_mat2_db = FT_PutMaterialInDB(mat = 'mat_2', comp_params = 'comp')
        
        # Put the parameters to build the interface in DB
        ft_interface_db = FT_PutInterfaceInDB(mat_1 = 'mat_1', mat_2 = 'mat_2',
                                             comp_params = 'comp',
                                             inter_params = 'inter')
        
        fw = Firework([ft_mat1, ft_mat2, ft_computation, ft_interface, 
                       ft_mat1_db, ft_mat2_db, ft_interface_db], name = fw_name)
        
        return fw

    @staticmethod
    def checkinp_homo_interface(material, computational, interface, 
                                fw_name='Check input parameters'):
        """
        Create a Fireworks to check if the necessary input for a homogeneous
        workflow are given and assignes default values to optional parameters. 
        Input parameters are checked for correct type and location in the spec.
        
        Parameters
        ----------
        material : list of str
            Keys list in fw_spec pointing to the input dict of the material.

        computational : list of str
            Keys list in fw_spec pointing to the computational parameters.

        interface : list of str
            Keys list in fw_spec pointing to the interface parameters, which
            are needed for creating the interface.

        fw_name : str
            Name of the returned FireWork
    
        Returns
        -------
        fw : fireworks.core.firework.Firework
            Firework that checks all input parameters for an heterogeneous WF.
            
        """
        
        # Firetasks checking the material parameters
        ft_mat = FTCheckInput(input_dict = material, read_key = 'material',
                               output_dict_name = 'mat')
        
        # Firetask checking the computational parameters
        ft_computation = FTCheckInput(input_dict = computational, 
                                      read_key = 'computational',
                                      output_dict_name = 'comp')
        
        # Firetask checking the interfacial matching parameters
        ft_interface = FTCheckInput(input_dict = interface,
                                    read_key = 'interface',
                                    output_dict_name = 'inter')
        
        # Put materials bulk and slab in DB
        ft_mat_db = FT_PutMaterialInDB(mat = 'mat', comp_params = 'comp')
        
        # Put the parameters to build the interface in DB
        ft_interface_db = FT_PutInterfaceInDB(mat_1 = 'mat', mat_2 = 'mat',
                                             comp_params = 'comp',
                                             inter_params = 'inter')
        
        fw = Firework([ft_mat, ft_computation, ft_interface, ft_mat_db, 
                       ft_interface_db], name = fw_name)
        
        return fw
