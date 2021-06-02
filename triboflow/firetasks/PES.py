#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:52:57 2020

@author: wolloch
"""
from operator import itemgetter
from monty.json import jsanitize

from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.core.operations import SymmOp
from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk

from triboflow.phys.high_symmetry import (
    get_slab_hs, get_interface_hs, pbc_hspoints, fix_hs_dicts)
from triboflow.phys.potential_energy_surface import get_pes
from triboflow.utils.plot_tools import plot_pes
from triboflow.utils.database import (
    Navigator, StructureNavigator, convert_image_to_bytes)
from triboflow.utils.vasp_tools import get_custom_vasp_relax_settings
from triboflow.utils.structure_manipulation import (
    interface_name, clean_up_site_properties, stack_aligned_slabs, 
    recenter_aligned_slabs)
from triboflow.workflows.base import dynamic_relax_swf

@explicit_serialize
class FT_StartPESCalcSubWF(FiretaskBase):
    """ Start a PES subworkflow.
    
    Starts a PES subworkflow using data from the high-level database.
    This is intended to be used to start a PES subworkflow from a main
    workflow.
    
    Parameters
    ----------
    mp_id_1 : str
        Materials Project database ID for the first material of the interface.
    mp_id_2 : str
        Materials Project database ID for the second material of the interface.
    miller_1 : list of int or str
        Miller index of the first material.
    miller_2 : list of int or str
        Miller index of the second material.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
        
    Returns
    -------
    FWAction that produces a detour PES subworkflow.    
    """
    required_params = ['mp_id_1', 'mp_id_2', 'miller_1', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'high_level_db']

    def run_task(self, fw_spec):

        from triboflow.workflows.subworkflows import calc_pes_swf
        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        miller_1 = self.get('miller_1')
        miller_2 = self.get('miller_2')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)
        
        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)
        
        nav_structure = StructureNavigator(
            db_file=db_file, 
            high_level=hl_db)
        interface_dict = nav_structure.get_interface_from_db(
            name=name,
            functional=functional
        )
        comp_params = interface_dict['comp_parameters']
        top_slab = Slab.from_dict(interface_dict['top_aligned'])
        bot_slab = Slab.from_dict(interface_dict['bottom_aligned'])
        already_done = interface_dict.get('relaxed_structure@min')
        
        if not already_done:
            SWF = calc_pes_swf(top_slab=top_slab,
                               bottom_slab=bot_slab,
                               interface_name=name,
                               functional=functional,
                               comp_parameters=comp_params,
                               output_dir=None)

            return FWAction(detours=SWF, update_spec=fw_spec)

        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_ComputePES(FiretaskBase):
    """ Compute the PES for a given interface, plot and save it.
    
    Uses the perviously computed energies for the unique high-symmetry points
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
        
        interpolation, E_list, pes_data, data, to_plot = get_pes(
            hs_all=all_hs,
            E=E_unique,
            cell=struct.lattice.matrix,
            to_fig=False)
        
        if file_output:
            data.dump('Computet_PES_data_'+name+'.dat')
            pes_data.dump('Interpolated_PES_data_'+name+'.dat')
            
        plot_pes(to_plot, cell*2, to_fig=name)
        plot_name = 'PES_' + str(name) + '.png'
        pes_image_bytes = convert_image_to_bytes('./'+plot_name)
        
        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional+'.interface_data',
            fltr={'name': name},
            new_values={'$set': {'PES.rbf': jsanitize(interpolation),
                                 'PES.all_energies': jsanitize(E_list),
                                 'PES.pes_data': jsanitize(pes_data),
                                 'PES.image': pes_image_bytes}})


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
        energy_list=[]
        calc_output = {}
        for s in lateral_shifts.keys():
            label = tag + '_' + s
            x_shift = lateral_shifts.get(s)[0][0]
            y_shift = lateral_shifts.get(s)[0][1]
            vasp_calc = nav.find_data(
                collection='tasks', 
                fltr={'task_label': label})
            energy = vasp_calc['output']['energy']
            struct = vasp_calc['output']['structure']
            energy_list.append([s, x_shift, y_shift, energy])
            calc_output[s] = {'energy': energy,
                              'relaxed_struct': struct,
                              'task_id': vasp_calc['_id']}

        sorted_energy_list = sorted(energy_list, key=itemgetter(3))
        
        min_stacking = sorted_energy_list[0][0]
        max_stacking = sorted_energy_list[-1][0]
        calc_min = nav.find_data(
            collection='tasks',
            fltr={'task_label': tag+'_'+min_stacking})
        calc_max = nav.find_data(
            collection='tasks',
            fltr={'task_label': tag+'_'+max_stacking})
        struct_min = calc_min['output']['structure']
        struct_max = calc_max['output']['structure']

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional+'.interface_data',
            fltr={'name': name},
            new_values={
                '$set': 
                    {'relaxed_structure@min': 
                         struct_min,
                     'relaxed_structure@max': 
                         struct_max,
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
        
        # Top slab needs to be mirrored to find the high symmetry points at the
        # interface.
        mirror = SymmOp.reflection(normal=[0,0,1], origin=[0, 0, 0])
        flipped_top = top_slab.copy()
        flipped_top.apply_operation(mirror)
        top_hsp_unique, top_hsp_all = get_slab_hs(flipped_top)
        
        bottom_hsp_unique, bottom_hsp_all = get_slab_hs(bot_slab)
        
        cell = bot_slab.lattice.matrix
        
        hsp_unique = get_interface_hs(bottom_hsp_unique, top_hsp_unique, cell)
        hsp_all = get_interface_hs(bottom_hsp_all, top_hsp_all, cell)
        
        c_hsp_u, c_hsp_a = fix_hs_dicts(hsp_unique, hsp_all,
                                      top_slab, bot_slab)
           
        b_hsp_u =  pbc_hspoints(bottom_hsp_unique, cell)
        b_hsp_a =  pbc_hspoints(bottom_hsp_all, cell)
        t_hsp_u =  pbc_hspoints(top_hsp_unique, cell)
        t_hsp_a =  pbc_hspoints(top_hsp_all, cell)
        
        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
            collection=functional+'.interface_data',
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
        
        inputs=[]
        for s in lateral_shifts.keys():
            label = tag + '_' + s
            x_shift = lateral_shifts.get(s)[0][0]
            y_shift = lateral_shifts.get(s)[0][1]
            #Make sure that there are no NoneTypes in the site_properties!
            inter_struct = stack_aligned_slabs(bot_slab,
                                             top_slab,
                                             top_shift = [x_shift, y_shift, 0])
            clean_struct = clean_up_site_properties(inter_struct)
            
            vis = get_custom_vasp_relax_settings(structure=clean_struct,
                                             comp_parameters=comp_params,
                                             relax_type='interface_z_relax')
            inputs.append([clean_struct, vis, label])
                        
        wf_name = 'PES relaxations for: '+name
        WF = dynamic_relax_swf(inputs_list = inputs, wf_name = wf_name)
        
        return FWAction(detours = WF)
