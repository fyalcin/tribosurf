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
from pymatgen.core.interface import Interface

from triboflow.phys.new_high_symm import InterfaceSymmetryAnalyzer
from triboflow.phys.new_potential_energy_surface import PESGenerator
from triboflow.phys.shaper import Shaper
from triboflow.utils.database import (
    Navigator, StructureNavigator, convert_image_to_bytes)
from triboflow.utils.plot_tools import plot_pes
from triboflow.utils.structure_manipulation import clean_up_site_properties
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
        all_shifts = inter_dict['PES']['high_symmetry_points']['all_shifts']
        unique_shifts = inter_dict['PES']['high_symmetry_points']['unique_shifts']
        energy_dict = inter_dict['PES']['high_symmetry_points']['energies_dict']
        group_assignments = inter_dict['PES']['high_symmetry_points']['group_assignments']
        
        if file_output:
            PG = PESGenerator(points_per_angstrom = 50,
                              interpolation_kernel = 'linear',
                              plot_hs_points = False,
                              plot_unit_cell = True,
                              plotting_ratio = 1.0,
                              normalize_minimum = True,
                              nr_of_contours = 30,
                              fig_title = name,
                              fig_type = 'png',
                              plot_path = './')
        else:
            PG = PESGenerator(points_per_angstrom = 50,
                              interpolation_kernel = 'linear',
                              plot_hs_points = False,
                              plot_unit_cell = True,
                              plotting_ratio = 1.0,
                              normalize_minimum = True,
                              nr_of_contours = 30,
                              fig_title = name,
                              fig_type = 'png',
                              plot_path = None)

        PG(interface=struct,
           energies_dict=energy_dict,
           all_shifts_dict=all_shifts,
           unique_shifts_dict=unique_shifts,
           group_names_dict=group_assignments)

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional + '.interface_data',
            fltr={'name': name},
            new_values={'$set': {'PES.rbf': jsanitize(PG.rbf),
                                 'PES.all_energies': jsanitize(PG.extended_energies),
                                 'PES.pes_data': jsanitize(PG.PES_on_meshgrid),
                                 'PES.image': PG.PES_as_bytes,
                                 'corrugation': PG.corrugation,
                                 'hsp@min': PG.hsp_min,
                                 'hsp@max': PG.hsp_max}})


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
        lateral_shifts = interface_dict['PES']['high_symmetry_points']['unique_shifts']
        group_assignments = interface_dict['PES']['high_symmetry_points']['group_assignments']

        nav = Navigator(db_file=db_file)
        ref_struct = Interface.from_dict(interface_dict['unrelaxed_structure'])
        area = np.linalg.norm(np.cross(ref_struct.lattice.matrix[0], ref_struct.lattice.matrix[1]))

        energy_list = []
        energy_dict = {}
        calc_output = {}
        for s in lateral_shifts.keys():
            label = tag + '_' + s
            vasp_calc = nav.find_data(
                collection='tasks',
                fltr={'task_label': label})
            struct = vasp_calc['output']['structure']
            energy = vasp_calc['output']['energy']
            energy *= 16.02176565 / area
            energy_list.append([s, group_assignments[s], energy])
            energy_dict[s] = energy
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
                     'PES.high_symmetry_points.energies_dict':
                         energy_dict}})


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
    interface : pymatgen.core.interface.Interface
        Interface object for which the PES is to be constructed
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

    required_params = ['interface', 'functional', 'interface_name']
    optional_params = ['db_file' 'high_level_db']

    def run_task(self, fw_spec):
        interface = self.get('interface')
        name = self.get('interface_name')
        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)


        ISA = InterfaceSymmetryAnalyzer()
        hsp_dict = ISA(interface)

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional + '.interface_data',
            fltr={'name': name},
            new_values={'$set':
                            {'PES.high_symmetry_points':
                                 {'bottom_unique': hsp_dict['bottom_high_symm_points_unique'],
                                  'bottom_all': hsp_dict['bottom_high_symm_points_all'],
                                  'top_unique': hsp_dict['top_high_symm_points_unique'],
                                  'top_all': hsp_dict['top_high_symm_points_all'],
                                  'unique_shifts': hsp_dict['unique_shifts'],
                                  'all_shifts': hsp_dict['all_shifts'],
                                  'group_assignments': hsp_dict['group_assignments']}}},
            upsert=True)

        return FWAction(update_spec=({'lateral_shifts': hsp_dict['unique_shifts']}))


@explicit_serialize
class FT_StartPESCalcs(FiretaskBase):
    """Start z-relaxations for different lateral positions of an interface.
    
    Take a list of lateral shifts from the fw_spec and start relaxations
    for each one of them as parallel detours.
    heterogeneous_wf
    Parameters
    ----------
    interface : pymatgen.core.interface.Interface
        Interface object for which the PES is to be constructed
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

    required_params = ['interface', 'interface_name',
                       'comp_parameters', 'tag']
    optional_params = ['db_file']

    def run_task(self, fw_spec):

        interface = self.get('interface')
        interface_name = self.get('interface_name')
        comp_params = self.get('comp_parameters')
        tag = self.get('tag')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        lateral_shifts = fw_spec.get('lateral_shifts')
        if not lateral_shifts:
            raise SystemExit('Lateral shifts not found in the fw_spec./n'
                             'Please check your Firework for errors!')

        inputs = []
        for name, shift in lateral_shifts.items():
            label = tag + '_' + name
            interface.in_plane_offset = shift
            clean_struct = clean_up_site_properties(interface)

            vis = get_custom_vasp_relax_settings(structure=clean_struct,
                                                 comp_parameters=comp_params,
                                                 relax_type='interface_z_relax')
            inputs.append([clean_struct, vis, label])

        wf_name = 'PES relaxations for: ' + interface_name
        WF = dynamic_relax_swf(inputs_list=inputs,
                               wf_name=wf_name,
                               add_static=True)

        return FWAction(detours=WF)
