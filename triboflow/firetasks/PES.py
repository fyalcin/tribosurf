#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:52:57 2020

@author: wolloch
"""
from operator import itemgetter

import numpy as np
from atomate.utils.utils import env_chk
from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from monty.json import jsanitize
from pymatgen.core.structure import Structure

from triboflow.phys.high_symmetry import (
    get_slab_hs, get_interface_hs, pbc_hspoints, fix_hs_dicts)
from triboflow.phys.interface_matcher import flip_slab
from triboflow.phys.potential_energy_surface import get_pes
from triboflow.phys.shaper import Shaper
from triboflow.utils.database import (
    Navigator, StructureNavigator, convert_image_to_bytes)
from triboflow.utils.plot_tools import plot_pes
from triboflow.utils.structure_manipulation import (
    clean_up_site_properties, stack_aligned_slabs,
    recenter_aligned_slabs)
from triboflow.utils.vasp_tools import get_custom_vasp_relax_settings
from triboflow.workflows.base import dynamic_relax_swf


@explicit_serialize
class FT_ComputePES(FiretaskBase):
    """ Compute the PES for a given interface, plot and save it.
    
    Uses the previously computed energies for the unique high-symmetry points
    and copies them to all the correct replica points. Replicates the points
    and fits the PES using radial basis functions. Output is saved in the
    database and if wanted also to files.
    
    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    file_output : bool
        Determines if results are written to disc.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.

    """

    required_params = ['interface_name', 'functional', 'file_output']
    optional_params = ['db_file', 'high_level_db']

    def run_task(self, fw_spec):

        name = self.get('interface_name')
        functional = self.get('functional')
        file_output = self.get('file_output')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        inter_dict = nav_structure.get_interface_from_db(
            name=name,
            functional=functional
        )
        struct = Structure.from_dict(inter_dict['relaxed_structure@min'])

        # Copy the energies for the unique points to all points
        E_unique = inter_dict['PES']['high_symmetry_points']['energy_list']
        all_hs = inter_dict['PES']['high_symmetry_points']['combined_all']
        cell = struct.lattice.matrix

        hsp_min = E_unique[0][0]
        hsp_max = E_unique[-1][0]

        interpolation, E_list, pes_data, data, to_plot = get_pes(
            hs_all=all_hs,
            E=E_unique,
            cell=struct.lattice.matrix,
            to_fig=False)

        corrugation = max(pes_data[:, 2]) - min(pes_data[:, 2])

        if file_output:
            data.dump('Computed_PES_data_' + name + '.dat')
            pes_data.dump('Interpolated_PES_data_' + name + '.dat')

        plot_pes(to_plot, cell, to_fig=name)
        plot_name = 'PES_' + str(name) + '.png'
        pes_image_bytes = convert_image_to_bytes('./' + plot_name)

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional + '.interface_data',
            fltr={'name': name},
            new_values={'$set': {'PES.rbf': jsanitize(interpolation),
                                 'PES.all_energies': jsanitize(E_list),
                                 'PES.pes_data': jsanitize(pes_data),
                                 'PES.image': pes_image_bytes,
                                 'corrugation': corrugation,
                                 'hsp@min': hsp_min,
                                 'hsp@max': hsp_max}})


@explicit_serialize
class FT_RetrievePESEnergies(FiretaskBase):
    """Retrieve the energies from the PES relaxations and update the db.
    
    Uses a tag together with the labels of the high-symmetry points saved
    in the high level database to retrieve the correct energies for each
    lateral shift of the interface. Sort the shifts by energies and save both
    the configuration with the lowest and the highest energy in the high level
    database. Also save the list of shifts and energies with corresponding
    labels there.

    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    tag : str
        Unique tag to identify the calculations.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
        
    Returns
    -------
    FWActions that produce a detour workflow with relaxations for the PES.
    """

    required_params = ['interface_name', 'functional', 'tag']
    optional_params = ['db_file', 'high_level_db']

    def run_task(self, fw_spec):

        name = self.get('interface_name')
        functional = self.get('functional')
        tag = self.get('tag')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        interface_dict = nav_structure.get_interface_from_db(
            name=name,
            functional=functional)
        lateral_shifts = interface_dict['PES']['high_symmetry_points']['combined_unique']

        nav = Navigator(db_file=db_file)
        ref_struct = Structure.from_dict(nav.find_data(
            collection='tasks',
            fltr={'task_label': f'{tag}_{list(lateral_shifts)[0]}'})['output']['structure'])
        area = np.linalg.norm(np.cross(ref_struct.lattice.matrix[0], ref_struct.lattice.matrix[1]))

        energy_list = []
        calc_output = {}
        for s in lateral_shifts.keys():
            label = tag + '_' + s
            x_shift = lateral_shifts.get(s)[0][0]
            y_shift = lateral_shifts.get(s)[0][1]
            vasp_calc = nav.find_data(
                collection='tasks',
                fltr={'task_label': label})
            struct = vasp_calc['output']['structure']
            energy = vasp_calc['output']['energy']
            energy *= 16.02176565 / area
            energy_list.append([s, x_shift, y_shift, energy])
            calc_output[s] = {'energy': energy,
                              'relaxed_struct': struct,
                              'task_id': vasp_calc['_id']}

        sorted_energy_list = sorted(energy_list, key=itemgetter(3))

        min_stacking = sorted_energy_list[0][0]
        max_stacking = sorted_energy_list[-1][0]
        calc_min = nav.find_data(
            collection='tasks',
            fltr={'task_label': tag + '_' + min_stacking})
        calc_max = nav.find_data(
            collection='tasks',
            fltr={'task_label': tag + '_' + max_stacking})
        struct_min_dict = calc_min['output']['structure']
        struct_max_dict = calc_max['output']['structure']

        struct_min = Structure.from_dict(calc_min['output']['structure'])
        struct_max = Structure.from_dict(calc_max['output']['structure'])
        inter_dist_min = Shaper.get_layer_spacings(struct_min)[0]
        inter_dist_max = Shaper.get_layer_spacings(struct_max)[0]

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional + '.interface_data',
            fltr={'name': name},
            new_values={
                '$set':
                    {'relaxed_structure@min':
                         struct_min_dict,
                     'relaxed_structure@max':
                         struct_max_dict,
                     'interface_distance@min':
                         inter_dist_min,
                     'interface_distance@max':
                         inter_dist_max,
                     'PES.calculations':
                         calc_output,
                     'PES.high_symmetry_points.energy_list':
                         sorted_energy_list}})


@explicit_serialize
class FT_FindHighSymmPoints(FiretaskBase):
    """Compute high symmetry points for the top and bottom slab and the interface.
    
    Finds the high symmetry points of the top side of the bottom slab and the
    bottom side of the top slab. This is done twice, once omitting duplicates,
    and once allowing them. It is made sure that the results are cartesian
    coordinates that lie inside the unit cell. The lists are combined so that
    every combination of unique points for the interface is present. The PES
    section of the high level database is updated with the results and the
    fw_spec is updated with the lateral shifts needed for the PES relaxations
    as well.
    
    Parameters
    ----------
    top_slab : pymatgen.core.surface.Slab
        Top slab of the interface.
    bottom_slab : pymatgen.core.surface.Slab
        Bottom slab of the interface.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    interface_name : str
        Name of the interface in the high-level database.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
        
    Returns
    -------
    FWActions that updates the fw_spec with lateral shifts.
    """

    required_params = ['top_slab', 'bot_slab', 'functional', 'interface_name']
    optional_params = ['db_file' 'high_level_db']

    def run_task(self, fw_spec):
        top_slab = self.get('top_slab')
        bot_slab = self.get('bot_slab')
        name = self.get('interface_name')
        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        # Top slab needs to be flipped to find the high symmetry points at the
        # interface.
        flipped_top = flip_slab(top_slab)
        top_hsp_unique, top_hsp_all = get_slab_hs(flipped_top)

        bottom_hsp_unique, bottom_hsp_all = get_slab_hs(bot_slab)

        cell = bot_slab.lattice.matrix

        hsp_unique = get_interface_hs(bottom_hsp_unique, top_hsp_unique, cell)
        hsp_all = get_interface_hs(bottom_hsp_all, top_hsp_all, cell)

        c_hsp_u, c_hsp_a = fix_hs_dicts(hsp_unique, hsp_all,
                                        top_slab, bot_slab)

        b_hsp_u = pbc_hspoints(bottom_hsp_unique, cell)
        b_hsp_a = pbc_hspoints(bottom_hsp_all, cell)
        t_hsp_u = pbc_hspoints(top_hsp_unique, cell)
        t_hsp_a = pbc_hspoints(top_hsp_all, cell)

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional + '.interface_data',
            fltr={'name': name},
            new_values={'$set':
                            {'PES.high_symmetry_points':
                                 {'bottom_unique': b_hsp_u,
                                  'bottom_all': b_hsp_a,
                                  'top_unique': t_hsp_u,
                                  'top_all': t_hsp_a,
                                  'combined_unique': jsanitize(c_hsp_u),
                                  'combined_all': jsanitize(c_hsp_a)}}},
            upsert=True)

        return FWAction(update_spec=({'lateral_shifts': c_hsp_u}))


@explicit_serialize
class FT_StartPESCalcs(FiretaskBase):
    """Start z-relaxations for different lateral positions of an interface.
    
    Take a list of lateral shifts from the fw_spec and start relaxations
    for each one of them as parallel detours.
    heterogeneous_wf
    Parameters
    ----------
    top_slab : pymatgen.core.surface.Slab
        Top slab of the interface.
    bottom_slab : pymatgen.core.surface.Slab
        Bottom slab of the interface.
    interface_name : str
        Name of the interface in the high-level database.
    comp_parameters : dict
        Computational parameters to be passed to the vasp input file generation.
    tag : str
        Unique tag to identify the calculations.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
        
    Returns
    -------
    FWActions that produce a detour workflow with relaxations for the PES.
    """

    required_params = ['top_slab', 'bot_slab', 'interface_name',
                       'comp_parameters', 'tag']
    optional_params = ['db_file']

    def run_task(self, fw_spec):

        top_slab = self.get('top_slab')
        bot_slab = self.get('bot_slab')
        name = self.get('interface_name')
        comp_params = self.get('comp_parameters')
        tag = self.get('tag')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        lateral_shifts = fw_spec.get('lateral_shifts')
        if not lateral_shifts:
            raise SystemExit('Lateral shifts not found in the fw_spec./n'
                             'Please check your Firework for errors!')

        top_slab, bot_slab = recenter_aligned_slabs(top_slab, bot_slab)

        # # List all sites of interface that have positive c coordinates as they
        # # are in the upper slab.
        # sites_to_shift = []
        # for i, s in enumerate(struct.sites):
        #     if s.c > 0:
        #         sites_to_shift.append(i)

        inputs = []
        for s in lateral_shifts.keys():
            label = tag + '_' + s
            x_shift = lateral_shifts.get(s)[0][0]
            y_shift = lateral_shifts.get(s)[0][1]
            # Make sure that there are no NoneTypes in the site_properties!
            inter_struct = stack_aligned_slabs(bot_slab,
                                               top_slab,
                                               top_shift=[x_shift, y_shift, 0])
            clean_struct = clean_up_site_properties(inter_struct)

            vis = get_custom_vasp_relax_settings(structure=clean_struct,
                                                 comp_parameters=comp_params,
                                                 relax_type='interface_z_relax')
            inputs.append([clean_struct, vis, label])

        wf_name = 'PES relaxations for: ' + name
        WF = dynamic_relax_swf(inputs_list=inputs,
                               wf_name=wf_name,
                               add_static=True)

        return FWAction(detours=WF)
