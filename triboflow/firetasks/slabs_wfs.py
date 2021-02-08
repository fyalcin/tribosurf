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
from uuid import uuid4
from pymatgen.core.structure import Structure
from fireworks import Workflow, Firework

from triboflow.utils.database import Navigator, NavigatorMP
from triboflow.firetasks.slabs_wfs import SlabWFs, SlabThicknessError
from triboflow.firetasks.slabs import FT_OptimalThickness
from triboflow.firetasks.errors import ReadSubWFsError
from triboflow.tasks.io import read_json

currentdir = os.path.dirname(__file__)

# ============================================================================
# Classes
# ============================================================================

class SlabWFs:
    """
    Collection of static methods to manipulate slab structures.
    
    """

    @staticmethod
    def conv_slabthick_surfene(structure, mp_id, miller, functional='PBE',
                               comp_params={}, spec={}, db_file=None,
                               low_level = None, high_level = 'triboflow', 
                               bulk_name = None, slab_name = None,
                               nplane_start=4, nplane_incr=2, nplane_conv=5, 
                               **kwargs):
        """
        Function to set the computational and physical parameters and start a 
        workflow to converge the thickness of the provided slabs.

        """                                                                 

        name = 'Slab Thickness optimization of ' + \
            structure.composition.reduced_formula + ' ' + str(miller)
        
        tag = name + '_' + str(uuid4())
        
        dfl = currentdir + '/defaults_fw.json'
        p = read_subwfs_params(default_file = dfl, default_key="subworkflow", **kwargs)

        # Check if the chemical formula passed is the same on MP database
        if mp_id.startswith('mp-') and flag[3:].isdigit():
            nav_mp = NavigatorMP()
            formula_from_struct = structure.composition.reduced_formula
            formula_from_flag = nav_mp.get_property_from_mp(flag, 
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
        if print_help:
            print('Once you workflow has finished you can access the '
                  'results from the database using this code:\n\n'
                  'import pprint\n'
                  'from triboflow.utils.database import GetBulkFromDB\n\n'
                  'nav = Navigator({})'
                  'results = find_data({} + ".slab_data", {"mpid": {}, "miller": {}})\n'
                  'pprint.pprint(results)\n'.format(db_file, functional, flag, miller))
        
        # Create a Workflow to calculate the Optimal Thickness for a structure
        ft_opt_thick = FT_OptimalThickness(structure,
                                           flag,
                                           miller,
                                           comp_parameters,
                                           spec,
                                           functional,
                                           layer_start,
                                           layer_incr, 
                                           n_converge, 
                                           db_file, 
                                           file_output,
                                           output_dir, 
                                           remote_copy,
                                           server,
                                           user,
                                           port,
                                           print_help)
        
        fw = Firework(ft_opt_thick, spec=spec, name='Converge Slab Thickness')
        wf = Workflow([fw], name=name)

        return wf

    @staticmethod
    def conv_slabthick_alat():        
        """
        Method description
        """
        pass


def read_subwfs_params(default_file, default_key, **kwargs):
    """
    [summary]

    Parameters
    ----------
    default_file : [type]
        [description]
    default_key : [type]
        [description]

    Returns
    -------
    [type]
        [description]

    Raises
    ------
    ReadSubWFsError
        [description]
    """

    defaults = read_json(default_file)
    defaults = defaults[default_key]
    params = {}

    if not set(kwargs.keys()).issubset(set(defaults.keys())):
        raise ReadSubWFsError("Values passed in kwargs are not known. Allowed "
                              "values: {}".format(defaults.keys()))

    # Set the parameters, passed by input or default values
    for key, value in defaults.items():
        params[key] = kwargs.get(key, value)
    
    return params