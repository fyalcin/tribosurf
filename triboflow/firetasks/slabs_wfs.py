#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:02:48 2021

Classes to generate subworkflows for calculating slab properties.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

#from fireworks import FiretaskBase, FWAction
#from atomate.utils.utils import env_chk

import os
from uuid import uuid4

from fireworks import Workflow, Firework

from triboflow.utils.database import Navigator, NavigatorMP
from triboflow.firetasks.slabs import (
    FT_StartThickConvo,
    FT_EndThickConvo
)
from triboflow.utils.errors import SlabThicknessError, SubWFError
from triboflow.tasks.io import read_json

currentdir = os.path.dirname(__file__)

# ============================================================================
# Classes
# ============================================================================

class SlabWF:
    """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

    Collection of static methods to manipulate slab structures.
    
    """

    @staticmethod
    def make_and_relax_slab(structure, mp_id, miller, functional,):
        pass

    @staticmethod
    def conv_slabthick_surfene(structure, mp_id, miller, functional='PBE',
                               comp_params={}, spec={}, db_file=None,
                               low_level=None, high_level='triboflow',
                               relax_type='slab_pos_relax', thick_min=4, 
                               thick_max=12, thick_incr=2, vacuum=10,
                               in_unit_planes=True, ext_index=0,
                               parallelization='low', recursion=False,
                               cluster_params={}):
        """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

        Function to set the computational and physical parameters and start a 
        workflow to converge the thickness of the provided slabs.

        """     

        # Set the Workflow name
        name = 'Slab Thickness optimization of ' + \
                structure.composition.reduced_formula + ' ' + str(miller)
        
        # Set the cluster parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_default_params(default_file=dfl, 
                                default_key="cluster_params", 
                                cluster_params=cluster_params)
        
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
                                                  convo_kind='surfene',
                                                  relax_type=relax_type,
                                                  thick_min=thick_min,
                                                  thick_max=thick_max,
                                                  thick_incr=thick_incr,
                                                  vacuum=vacuum,
                                                  in_unit_planes=in_unit_planes,
                                                  ext_index=ext_index,
                                                  parallelization=parallelization,
                                                  recursion=recursion,
                                                  cluster_params=p)

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
                                              parallelization=parallelization, 
                                              recursion=recursion,
                                              cluster_params=cluster_params)

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
        """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

        Subworkflow to converge the slab thickness by converging the lattice
        parameters. Not implemented yet.
        """
        pass
    
    @staticmethod
    def _check_subwf_params(self, structure, mp_id, miller, functional, db_file, 
                            comp_params, cluster_params):
        """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

        Check if the parameters passed to the Firetasks are correct or not and
        print information.
        ** This is a temporary method, to be refactored more logically. **
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
        if cluster_params['print_help']:
            print('Once you workflow has finished you can access the '
                  'results from the database using this code:\n\n'
                  'import pprint\n'
                  'from triboflow.utils.database import GetBulkFromDB\n\n'
                  'nav = Navigator({})'
                  'results = find_data({} + ".slab_data", dict("mpid": {}, '
                  '"miller": {}))\n'
                  'pprint.pprint(results)\n'.format(db_file, functional, mp_id, miller))

def read_default_params(default_file, default_key, dict_params):

    defaults = read_json(default_file)
    defaults = defaults[default_key]
    params = {}

    if not set(dict_params.keys()).issubset(set(defaults.keys())):
        raise SubWFError("The values passed as dictionary params are not known. "
                         "Allowed values for {}, read in {}: {}".format(
                             default_key, default_file, defaults.keys()))

    # Set the parameters, passed by input or default values
    for key, value in defaults.items():
        params[key] = dict_params.get(key, value)
    
    return params
