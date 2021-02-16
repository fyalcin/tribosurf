#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:47:15 2021

Workflows for the surface energy

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 9th, 2021'

from uuid import uuid4

import numpy as np
from pymatgen.core.surface import SlabGenerator
from fireworks import Workflow, Firework

from triboflow.firetasks.slabs import FT_GenerateSlabs, FT_RelaxStructure
from triboflow.utils.errors import SubWFError


class SurfEneWF:

    @staticmethod
    def surface_energy(structure, mp_id, miller, functional='PBE', spec={},
                       db_file=None, database=None, relax_type="slab_pos_relax",
                       comp_params={}, thick_min=4, thick_max=12, thick_incr=2, 
                       vacuum=10, in_unit_planes=True, ext_index=0, 
                       cluster_params={}, parallelization='low', recursion=False):
        """
        Description of the method...

        """

        # Define the first Firework
        # =====================================================================

        # Define the thicknesses to create and simulate the slabs
        if parallelization == 'low':
            if recursion:
                thickness = [thick_min]
            else:
                thickness = [0, thick_min, thick_max]

        elif parallelization == 'high':
            try:
                thickness = np.arange(thick_min, thick_max, thick_incr)
                thickness = thickness.tolist()
                assert len(thickness)>0

            except: SubWFError("Thickness parameters are not defined correctly")
        
        else:
            raise SubWFError("The value passed as parallelization is not known "
                              "Allowed value: 'low', 'high'")

        # Create the dictionary key where the unrelaxed slab will be saved
        miller_str = get_miller_str(miller)
        # formula = structure.composition.reduced_formula
        slab_name = [['slab_' + miller_str + '_' + str(t), 'unrelaxed'] for t in thickness]

        # Generate the slabs and store them in the low level database, under
        # the dictionary key: 'unrelaxed_slab'
        ft_gen_slabs = FT_GenerateSlabs(structure=structure,
                                        mp_id=mp_id,
                                        miller=miller,
                                        collection=functional+'.slab_data',
                                        db_file=db_file,
                                        database=database,
                                        thickness=thickness,
                                        thick_max=thick_max,
                                        vacuum=vacuum,
                                        ext_index=ext_index,
                                        in_unit_planes=in_unit_planes,
                                        slab_name=slab_name)

        fw_gen_slabs = Firework([ft_gen_slabs],
                                spec = spec,
                                name = 'Generate ' + str(len(thickness)) + ' slabs')

        # Define the second Firework
        # ==================================================

        # Create the tags to store the calculation of the slabs
        tags = create_tags(slab_name)

        # Create the Firetasks to relax the structures
        fw_relax_slabs = []
        for i, name in enumerate(slab_name):
            ft = FT_RelaxStructure(mp_id=mp_id,
                                   functional=functional,
                                   collection=functional+'.slab_data',
                                   name=slab_name[i][0],
                                   db_file=db_file,
                                   database=database,
                                   relax_type=relax_type,
                                   comp_params=comp_params,
                                   miller=miller,
                                   check_key='relaxed',
                                   tag=tags[i])

            fw = Firework([ft],
                          spec = spec,
                          name = 'Relax slab: ' + name)
            fw_relax_slabs.append(fw)

        # Define the third Firework(s)
        # ==================================================
        ft_surfene = FT_SurfaceEnergy()

        # Define the name of the slab(s)


        fw = Firework([ft1, ft2],
                        spec = spec,
                        name = 'Calculate the Surface Energy')
        
        Workflow()

def create_tags(name):

    # Create a list of tags
    if isinstance(name, list):
        tag = [n + '_' + str(uuid4()) for n in name]

    else:
        tag = name + '_' + str(uuid4())
    
    return tag

def get_miller_str(miller):
    return ''.join(str(s) for s in miller)