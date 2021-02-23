#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:47:15 2021

Workflows for calculating and converging the surface energy of a slab.

The module contains:
    
** SurfEneWF **:
    General class to work on crystalline slabs, workflows are static method.
    It includes the following methods:
        - conv_slabthick_surfene
        - conv_slabthick_alat
        - _check_subwf_params
        
    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 9th, 2021'


import numpy as np
from fireworks import Workflow, Firework

from triboflow.utils.utils import create_tags, get_miller_str
from triboflow.firetasks.slabs import (
    FT_GenerateSlabs, 
    FT_RelaxStructure, 
    FT_MoveTagResults
)
from triboflow.firetasks.surfene import FT_SurfaceEnergy
from triboflow.utils.errors import SubWFError


class SurfEneWF:
    """
    Collection of static methods to calculate the surface energy of a slab. 
    It also contains methods to converge the surface energy of a slab.
    
    """

    @staticmethod
    def conv_surface_energy(structure, mp_id, miller, functional='PBE', spec={},
                            db_file=None, low_level=None, high_level='triboflow',
                            relax_type='slab_pos_relax', comp_params={}, thick_min=4, 
                            thick_max=12, thick_incr=2, vacuum=10, in_unit_planes=True, 
                            ext_index=0, cluster_params={}, parallelization='low', 
                            recursion=False):
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
        formula = structure.composition.reduced_formula
        miller_str = get_miller_str(miller)
        slab_entry = [[formula + 'slab_' + miller_str + '_' + str(t), 
                       'unrelaxed'] for t in thickness]

        # Generate the slabs and store them in the low level database, under
        # the dictionary key: 'unrelaxed_slab'
        ft_gen_slabs = FT_GenerateSlabs(structure=structure,
                                        mp_id=mp_id,
                                        miller=miller,
                                        collection=functional+'.slab_data',
                                        db_file=db_file,
                                        database=low_level,
                                        thickness=thickness,
                                        thick_bulk=thick_max,
                                        vacuum=vacuum,
                                        ext_index=ext_index,
                                        in_unit_planes=in_unit_planes,
                                        entry=slab_entry)

        fw_gen_slabs = Firework([ft_gen_slabs],
                                spec=spec,
                                name='Generate ' + str(len(thickness)) + ' slabs')

        # Define the second Firework
        # ==================================================

        # Create the tags to store the calculation of the slabs
        tags = create_tags(slab_entry)

        # Create the Firetasks to relax the structures
        fw_relax_slabs = []
        for thk, n, t in zip(thickness, slab_entry, tags):
            ft_1 = FT_RelaxStructure(mp_id=mp_id,
                                     functional=functional,
                                     collection=functional+'.slab_data',
                                     entry=n,
                                     tag=t,
                                     db_file=db_file,
                                     database=low_level,
                                     relax_type=relax_type,
                                     comp_params=comp_params,
                                     miller=miller,
                                     check_key='relaxed')

            ft_2 = FT_MoveTagResults(mp_id=mp_id,
                                     collection_from=functional+'.slab_data',
                                     collection_to=functional+'.slab_data',
                                     db_file=db_file,
                                     database_from=low_level,
                                     database_to=high_level,
                                     miller=miller,
                                     entry_check=[
                                         ['thickness', 
                                         'data_' + str(thk), 
                                         'calc_output']
                                         ],
                                     entry_to=[
                                         ['thickness', 
                                          'data_' + str(thk), 
                                          'calc_output'] * 9
                                         ],
                                     entry_from=[
                                         ['output', 'structure'],
                                         ['nsites'],
                                         ['output', 'density'],
                                         ['output', 'energy'],
                                         ['output', 'energy_per_atom' ],
                                         ['output', 'bandgap'],
                                         ['output', 'forces'],
                                         ['output', 'stresses'],                                         
                                         ['_id']
                                         ],
                                     struct_kind='slab',
                                     override=True,
                                     cluster_params=cluster_params)

            fw = Firework([ft_1, ft_2],
                          spec=spec,
                          name='Relax and store in DB, slab: ' + n)
            fw_relax_slabs.append(fw)

        # Define the third Firework(s)
        # ==================================================

        # Set the location of the energies in the high level DB
        entry_surfene = [['thickness', 'data_0', 'calc_output']]
        thk_loop = thickness if recursion else thickness[1:]
        for thk in thk_loop:
            entry_surfene.append(['thickness', 'data_' + str(thk), 'calc_output'])

        ft_surfene = FT_SurfaceEnergy(mp_id=mp_id,
                                      collection=functional+'.slab_data',
                                      miller=miller,
                                      entry=entry_surfene,
                                      db_file=db_file,
                                      database=high_level)

        fw_surfene = Firework([ft_surfene],
                              spec=spec,
                              name='Calculate the Surface Energies')

        # Define and return the Workflow
        # ==================================================       

        # Build the workflow list and the links between elements
        wf_list = [fw_gen_slabs, ]
        links = {fw_gen_slabs: fw_relax_slabs}
        for fw in fw_relax_slabs:
            wf_list.append(fw)
            links.update({fw: [fw_surfene]})
        wf_list.append(fw_surfene)

        wf = Workflow(wf_list, links, name="Converge surface energy WF")

        return wf
