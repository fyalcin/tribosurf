#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:02:48 2021

Collection of Workflow to study the properties and structures of crystal slabs.

The module contains:

** SlabWF **:
    General class to work on crystalline slabs, workflows are static method.
    It includes the following methods:
        - conv_slabthick_surfene
        - conv_slabthick_alat
        - _check_subwf_params

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'


import os

from fireworks import Workflow, Firework

from triboflow.utils.database import NavigatorMP
from triboflow.firetasks.run_slabs_wfs import (
    FT_StartThickConvo,
    FT_EndThickConvo
)
from triboflow.utils.utils import read_default_params


currentdir = os.path.dirname(__file__)
currentdir += '/../firetasks'


# ============================================================================
# Classes
# ============================================================================

class SlabWF:
    """
    Collection of static methods to create, cut, modify, merge, and generally
    manipulate slab structures. It also contains methods to manage the 
    optimization of a slab thickness, based on the evaluation of the surface 
    energy and the lattice parameter.

    """
    
    @staticmethod
    def conv_slabthick_surfene(structure, mp_id, miller, functional='PBE',
                               comp_params={}, spec={}, db_file=None,
                               low_level=None, high_level='triboflow',
                               relax_type='slab_pos_relax', thick_min=4, 
                               thick_max=12, thick_incr=2, vacuum=10,
                               in_unit_planes=True, ext_index=0, conv_thr=0.025,
                               parallelization='low', recursion=False,
                               cluster_params={}, override=False):
        """ 
        Function to set the computational and physical parameters and start a 
        workflow to converge the thickness of the provided slabs.

        """     

        # Set the name of the Workflow
        name = 'Slab Thickness optimization of ' + \
                structure.composition.reduced_formula + ' ' + str(miller)
        
        # Set the cluster parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_default_params(default_file=dfl, 
                                default_key="cluster_params", 
                                dict_params=cluster_params)
        
        # Print relevant information and raise errors based on parameters
        SlabWF._check_subwf_params(structure, mp_id, miller, functional, 
                                   db_file, comp_params, cluster_params)

        # Create a Firework to calculate the Optimal Thickness for a structure
        ft_start_thick_convo = FT_StartThickConvo(structure=structure,
                                                  mp_id=mp_id,
                                                  miller=miller,
                                                  functional=functional,
                                                  comp_params=comp_params,
                                                  db_file=db_file,
                                                  low_level=low_level,
                                                  high_level=high_level,
                                                  conv_kind='surfene',
                                                  relax_type=relax_type,
                                                  thick_min=thick_min,
                                                  thick_max=thick_max,
                                                  thick_incr=thick_incr,
                                                  vacuum=vacuum,
                                                  in_unit_planes=in_unit_planes,
                                                  ext_index=ext_index,
                                                  parallelization=parallelization,
                                                  recursion=recursion,
                                                  cluster_params=p,
                                                  override=override)

        ft_end_thick_convo = FT_EndThickConvo(structure=structure,
                                              mp_id=mp_id, 
                                              miller=miller, 
                                              functional=functional,
                                              comp_params=comp_params, 
                                              spec=spec,
                                              db_file=db_file,
                                              low_level=low_level, 
                                              high_level=high_level,
                                              relax_type=relax_type, 
                                              thick_min=thick_min,
                                              thick_max=thick_max, 
                                              thick_incr=thick_incr, 
                                              vacuum=vacuum,
                                              in_unit_planes=in_unit_planes, 
                                              ext_index=ext_index,
                                              conv_thr=conv_thr,
                                              parallelization=parallelization, 
                                              recursion=recursion,
                                              cluster_params=cluster_params,
                                              override=override)

        # Set it to a firework and a workflow
        # TODO: Understand if it is possible to have a structure of this kind
        # if a detour is done.
        fw = Firework([ft_start_thick_convo, ft_end_thick_convo],
                      spec = spec,
                      name = 'Converge slab thickness via surfene WF')           
        wf = Workflow([fw], name=name)

        return wf

    @staticmethod
    def conv_slabthick_alat():        
        """
        Subworkflow to converge the slab thickness by converging the lattice
        parameters. Not implemented yet.
        """
        pass

    @staticmethod
    def _check_subwf_params(structure, mp_id, miller, functional, db_file, 
                            comp_params, cluster_params):
        """
        Check if the Firetasks parameters are correct and print information.

        """

        # Check if the chemical formula passed is the same on MP database
        if mp_id.startswith('mp-') and mp_id[3:].isdigit():
            nav_mp = NavigatorMP()
            formula_from_struct = structure.composition.reduced_formula
            formula_from_flag = nav_mp.get_property_from_mp(mp_id, 
                                                            ['pretty_formula'])
            if not formula_from_flag == formula_from_struct:
                raise SystemExit('The chemical formula of your structure ({}) '
                                 'does not match the chemical formula of the flag '
                                 '(mp-id) you have chosen which corresponds '
                                 'to {}.\n'.format(
                                 formula_from_struct, formula_from_flag))

        # Check computational parameters and use defaults if necessary
        if comp_params == {}:
            print('\nNo computational parameters have been defined!\n'
                'Workflow will run with:\n'
                '   ISPIN = 1\n'
                '   ISMEAR = 0\n'
                '   ENCUT = 520\n'
                '   kpoint density kappa = 5000\n'
                'We recommend to pass a comp_parameters dictionary'
                ' of the form:\n'
                '   {"use_vdw": <True/False>,\n'
                '    "use_spin": <True/False>,\n'
                '    "is_metal": <True/False>,\n'
                '    "encut": <int>,\n'
                '    "k_dens": <int>}\n')

        # Print help to the user
        if'print_help' in cluster_params.keys():
            if cluster_params['print_help']:
                print('Once you workflow has finished you can access the '
                      'results from the database using this code:\n\n'
                      'import pprint\n'
                      'from triboflow.utils.database import GetBulkFromDB\n\n'
                      'nav = Navigator({})'
                      'results = find_data({} + ".slab_data", dict("mpid": {}, '
                      '"miller": {}))\n'
                      'pprint.pprint(results)\n'.format(db_file, functional, mp_id, miller))
