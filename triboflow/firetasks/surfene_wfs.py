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
    def surface_energy(structure, mp_id, miller, functional='PBE', 
                       comp_params={}, db_file=None, collection=None, 
                       thick_start=4, thick_incr=1, nsteps=6, vacuum=10, 
                       relax_type="slab_pos_relax", ext_index=True, 
                       slab_name=None, spec={}, cluster_params={}):
        """
        Description...

        Parameters
        ----------
        structure : [type]
            [description]
        mp_id : [type]
            [description]
        miller : [type]
            [description] ft_generate_slabs = FT_GenerateSlabs(structure=structure,
                                             mp_id=mp_id,
                                             miller=miller,
                                             functional=functional,
                                             db_file=db_file,
                                             collection=collection,
                                             thickness=thickness,
                                             vacuum=vacuum,
                                             ext_index=ext_index,
                                             symmetrize=False,
                                             slab_name
        functional : str, optional
            [description], by default 'PBE'
        db_file : [type], optional
            [description], by default None
        collection : [type], optional
            [description], by default None
        thick_start : int, optional
            [description], by default 4
        thick_incr : int, optional
            [description], by default 1
        nsteps : int, optional
            [description], by default 6
        vacuum : int, optional
            [description], by default 10
        ext_index : bool, optional
            [description], by default True
        slab_name : [type], optional
            [description], by default None
        spec : dict, optional
            [description], by default {}
        """

        # Define the first Firework
        # =====================================================================

        # Define the thicks that are required to calculate the surface energies
        thickness = np.arange(thick_start, thick_start + nsteps * thick_incr)
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
                                             collection=collection,
                                             thickness=thickness,
                                             vacuum=vacuum,
                                             ext_index=ext_index,
                                             symmetrize=False,
                                             slab_name=names)
        fw_genslabs = Firework([ft_generate_slabs],
                               spec = spec,
                               name = 'Generate ' + str(nsteps) + ' slabs' )

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
                                    collection=collection,
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
