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

import numpy as np
from fireworks import Workflow, Firework

from triboflow.firetasks.slabs import FT_GenerateSlabs, FT_RelaxStructure, FT_PutStructInDB


class SurfEneWF:

    @staticmethod
    def surface_energy(structure, mp_id, miller, functional='PBE', spec={},
                       db_file=None, database=None, comp_params={},  thick_min=4, 
                       thick_max=1, thick_incr=6, vacuum=10, cluster_params={},
                       relax_type="slab_pos_relax", ext_index=True, slab_name=None):
        """
        Description of the method...

        """

        # Define the first Firework
        # =====================================================================

        # Define the thicks that are required to calculate the surface energies
        thickness = np.arange(thick_min, thick_max, thick_incr)
        thickness = thickness.tolist()

        # Define the name of the slab(s)
        out_names = structure.composition.reduced_formula
        miller_str = ''.join(str(s) for s in miller)
        if slab_name is not None:
            out_names = slab_name
        names = [out_names + miller_str + '_' + str(t) for t in thickness]

        ft_generate_slabs = FT_GenerateSlabs(structure=structure,
                                             mp_id=mp_id,
                                             miller=miller,
                                             functional=functional,
                                             db_file=db_file,
                                             database=database,
                                             thickness=thickness,
                                             vacuum=vacuum,
                                             ext_index=ext_index,
                                             symmetrize=False,
                                             slab_name=names)
        fw_genslabs = Firework([ft_generate_slabs],
                               spec = spec,
                               name = 'Generate ' + str(len(thickness)) + ' slabs' )

        # Define the second Firework(s)
        # ==================================================
        fw_relax_slabs = []
        for slab_name in names:

            ft1 = FT_RelaxStructure(mp_id=mp_id,
                                    miller=miller,
                                    functional=functional,
                                    comp_params=comp_params,
                                    struct_kind='slab',
                                    db_file=db_file,
                                    database=database,
                                    name=slab_name,
                                    relax_type=relax_type)

            fw = Firework([ft1, ft2],
                          spec = spec,
                          name = 'Relax slab: ' + slab_name)

        # Define the third Firework(s)
        # ==================================================
        ft_surfene = FT_SurfaceEnergy()

        fw = Firework([ft1, ft2],
                        spec = spec,
                        name = 'Calculate the Surface Energy')
        
        Workflow()
