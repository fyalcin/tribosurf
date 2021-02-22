#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:16:55 2021

Collection of Firetasks to start workflows concerning slabs structure.
They can be used to run existing workflows as subworkflows within other ones.

The module contains the following Firetasks:

** Slab thickness convergence **

- FT_SlabOptThick 
    Starts a subworkflow to perform a convergence process to find the optimal 
    thickness of a slab structure. First step to be called within a workflow.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

- FT_StartThickConvo
    Start a subworkflow to select a desired convergence criterion to calculate
    the slab thickness. Implemented criteria: surface energy.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

- FT_EndThickConvo
    Check the results and call recursively the slab thickness workflow it the
    convergence has not been achieved yet.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 22nd, 2021'


import os

import numpy as np
from pymatgen.core.surface import Structure
from fireworks import explicit_serialize, FiretaskBase, FWAction

from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.utils import (
    read_runtask_params,
    get_one_info_from_dict
)
from triboflow.utils.errors import SlabOptThickError


currentdir = os.path.dirname(__file__)


@explicit_serialize
class FT_SlabOptThick(FiretaskBase):
    """
    Start a subworkflow as a detour to calculate the optimal thickness for a 
    slab generated from a supplied bulk structure and a certain orientation 
    (miller index). The thickness can be converted either by evaluating how the 
    surface energy or the lattice parameters (not implemented) change with the 
    number of atomic planes.

    Parameters
    ----------
    mp_id : str
        MP-id of the structure from the MP database.
    
    miller : list of int, or str
        Miller indexes (h, k, l) to select the slab orientation.
        
    functional : str
        Functional for the pseudopotential to be adopted.

    db_file : str or None
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.
        
    low_level : str or None, optional
        Name of the table of the "low level" database, saved in db_file. The 
        intermediate calculations and raw data during will be saved here. If 
        nothing is passed the Firework database is used. The default is None.

    high_level : str, optional
        Name of the table of the "high level" database, saved in db_file.
        The slab optimal thickness will be saved here. The slab energy and
        surface energy will be saved too if conv_kind is 'surfene'.
        The default is 'triboflow'.

    conv_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene', 
        'alat'. The latter is not implemented yet. The default is 'surfene'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.

    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point. In
        case of low parallelization this gives the value of atomic layers for
        the only slab which is created and relaxed at each step after the first.
        The default is 4.
        
    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If convergence is
        not reached this value is considered the optimal one. The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step during
        the iterative procedure. The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.
    
    ext_index : int, optional
        Use the ext_index element from SlabGenerator.get_slabs as a slab.
        The default is 0.
    
    conv_thr : float, optional
        Threshold for the convergence process. The default is 0.025.

    bulk_entry : str or list or None, optional
        Name of the custom bulk dictionary, to be retrieved from the high level
        database and to be used to build the slabs. Bulks are identified by 
        mp_id and functional but there might be different structures of the
        same material. The default is "structure_fromMP".

    slab_entry : str or None, optional
        Where to search for the information about the optimal thickness in the
        slab dictionary in the high level database. The default is None.

    """
    
    _fw_name = 'Start a subworkflow to converge slab thickness'

    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'low_level', 'high_level', 'conv_kind',
                       'relax_type', 'thick_min', 'thick_max', 'thick_incr',
                       'vacuum', 'in_unit_planes', 'ext_index', 'conv_thr',
                       'bulk_entry', 'slab_entry']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, self.optional_params,
                                default_file=dfl, default_key="SlabOptThick")
        
        # Retrieve the bulk information from the high level DB
        nav_struct = StructureNavigator(p['db_file'], p['high_level'])
        slab = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                           functional=p['functional'], 
                                           miller=p['miller'])
        
        # If data is saved elsewhere from standard position
        if p['slab_entry'] is not None:
            slab = get_one_info_from_dict(slab, p['slab_entry'])

        # Start a subworkflow to converge the thickness if not already done
        stop_convergence = slab.get('opt_thickness', None)

        if not stop_convergence:
            
            # Retrieve the desired bulk structure
            bulk = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                               functional=p['functional'])
            bulk = get_one_info_from_dict(bulk, p['bulk_entry'])

            structure = Structure.from_dict(bulk.get('structure_fromMP'))
            comp_params = slab.get('comp_params', {})

            wf = self.select_slabthick_conv(structure=structure, 
                                            comp_params=comp_params,
                                            p=p)

            return FWAction(detours=wf, update_spec=fw_spec)

        # Continue the Workflow
        else:
            return FWAction(update_spec=fw_spec)

    def select_slabthick_conv(self, structure, comp_params, p):
        """
        Select the desired subworkflow from the SlabWFs class, to converge the 
        slab thickness either by evaluating the surface energy or the lattice 
        parameter.
        
        """
        
        from triboflow.workflows.slabs_wfs import SlabWF

        if p['conv_kind'] == 'surfene':
            generate_wf = SlabWF.conv_slabthick_surfene
        elif p['conv_kind'] == 'alat':
            generate_wf = SlabWF.conv_slabthick_alat
        else:
            raise SlabOptThickError("Wrong input argument for conv_kind. "
                                    "Allowed options: 'surfene', 'alat'")

        wf = generate_wf(structure=structure, mp_id=p['mp_id'], 
                         miller=p['miller'], functional=p['functional'], 
                         comp_params=p['comp_params'], db_file=p['db_file'],
                         low_level=p['low_level'], high_level=p['high_level'],
                         relax_type=p['relax_type'], thick_min=p['thick_min'], 
                         thick_max=p['thick_max'], thick_incr=p['thick_incr'], 
                         vacuum=p['vacuum'], in_unit_planes=p['in_unit_planes'],
                         ext_index=p['ext_index'], conv_thr=p['conv_thr'])

        return wf

@explicit_serialize
class FT_StartThickConvo(FiretaskBase):
    """
    It starts a subworkflow as a detour to converge the thickness of a slab.
    The thickness can be converged either by evaluating how the surface energy
    or the lattice parameters changes with the number of atomic planes.
    It is the first element of the SlabWF.conv_slabthick_surfene workflow.

    (At the moment only the surface energy convergence is implemented)
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Bulk structure to build the slabs.

    mp_id : str
        MP-id of the structure from the MP database.
    
    miller : list of int, or str
        Miller indexes (h, k, l) to select the slab orientation.
        
    functional : str
        Functional for the pseudopotential to be adopted.

    db_file : str or None
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.
        
    low_level : str or None, optional
        Name of the table of the "low level" database, saved in db_file. The 
        intermediate calculations and raw data during will be saved here. If 
        nothing is passed the Firework database is used. The default is None.

    high_level : str, optional
        Name of the table of the "high level" database, saved in db_file.
        The slab optimal thickness will be saved here. The slab energy and
        surface energy will be saved too if conv_kind is 'surfene'.
        The default is 'triboflow'.

    conv_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene', 
        'alat'. The latter is not implemented yet. The default is 'surfene'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.

    comp_params : dict, optional
        Computational parameters for the VASP simulations. If not set, default 
        parameters will be used instead. The default is {}.
    
    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point. In
        case of low parallelization this gives the value of atomic layers for
        the only slab which is created and relaxed at each step after the first.
        The default is 4.
        
    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If convergence is
        not reached this value is considered the optimal one. The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step during
        the iterative procedure. The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.

    ext_index : int, optional
        Use the ext_index element from SlabGenerator.get_slabs as a slab.
        The default is 0.

    cluster_params : dict, optional
        Dictionary containing cluster-related options to run efficiently the
        VASP simulations on a cluster. The default is {}.

    """
    
    _fw_name = 'Start the slab thickness convergence'
    required_params = ['structure', 'mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'low_level', 'high_level', 'conv_kind', 
                       'relax_type', 'comp_params', 'thick_min', 'thick_max', 
                       'thick_incr', 'vacuum', 'in_unit_planes', 'ext_index', 
                       'parallelization', 'recursion', 'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, 
                                self.optional_params, default_file=dfl, 
                                default_key="StartThickConvo")
        
        # Select the convergence of interest
        wf = self.select_conv(p)

        return FWAction(detours=wf, update_spec=fw_spec)

        def select_conv(self, p):
        
            from triboflow.workflows.surfene_wfs import SurfEneWF

            if p['conv_kind'] == 'surfene':
                wf = SurfEneWF.surface_energy(structure=p['structure'],
                                            mp_id=p['mp_id'], 
                                            miller=p['miller'], 
                                            functional=p['functional'],
                                            db_file=p['db_file'], 
                                            low_level=p['low_level'],
                                            high_level=p['high_level'],
                                            relax_type=p['relax_type'],
                                            comp_params=p['comp_params'],
                                            thick_min=p['thick_min'], 
                                            thick_max=p['thick_max'],
                                            thick_incr=p['thick_incr'],
                                            vacuum=p['vacuum'],
                                            in_unit_planes=p['in_unit_planes'],
                                            ext_index=p['ext_index'], 
                                            parallelization=p['parallelization'],
                                            recursion=p['recursion'],
                                            cluster_params=p['cluster_params'])
                return wf

            else:
                raise SystemExit('Lattice parameter convergence not yet implemented')

@explicit_serialize
class FT_EndThickConvo(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['structure', 'mp_id', 'miller']
    optional_params = ['db_file', 'low_level', 'high_level', 'functional', 
                       'conv_kind', 'relax_type', 'comp_params', 'thick_min', 
                       'thick_max', 'thick_incr', 'vacuum', 'in_unit_planes', 
                       'ext_index', 'conv_thr', 'parallelization', 'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, 
                                self.optional_params, default_file=dfl, 
                                default_key="EndThickConvo")
        case = self._get_case()
        
        # Retrieve the surface energy
        data, index = self.get_data(case, p)

        # Analyze the data, if convergence is reached data is stored in DB
        stop_convergence = self.analyze_data(data, index, case, p)

        # Decide whether to rerun recursively the surface energy workflow
        if not stop_convergence:
            wf = self.call_recursion()

            return FWAction(detours=wf, update_spec=fw_spec)

        else:
            return FWAction(update_spec=fw_spec)

        def _get_case(self, p):
            """
            
            """
            # Select the keys to be used when reading the dictionary
            if p['conv_kind'] == 'surfene':
                case = 'surface_energy'
            elif p['conv_kind'] == 'alat':
                case = 'lattice'
            else:
                raise SlabOptThickError("Wrong argument: 'conv_kind'. Allowed "
                                        "values: 'surfene', 'alat'")
            return case

        def get_data(self, case, p):
            """
            Extract the surface energies from the high level DB.

            """

            # Call the navigator for retrieving the thickness dict out of DB
            nav = Navigator(db_file=p['db_file'], high_level=p['high_level'])
            thickness_dict = nav.find_data(p['functional'] + '.slab_data', 
                                           {'mpid': p['mp_id'], 
                                           'miller': p['miller']})['thickness']
            
            # Get the indexes with the thickness and the desired data
            index = []
            data = []
            for key, item in thickness_dict.values():
                if key.startswith('data_'):
                    index.append(key.split('_')[-1])
                    data.append(item['calc_output', case])
            sorted_index = index.argsort()
            index = np.array(index[sorted_index])
            data = np.array(data[sorted_index])

            # Check for consistency with the value from the max thickness
            dmax = thickness_dict['data_' + str(p['max_thick'])]['calc_output'][case]
            if dmax != data[-1]:
                raise SlabOptThickError("An unexpected error occurred")
            
            return data, index
        
        def analyze_data(self, data, index, case, p):

            # Calculate the relative error to the last element
            error_to_last = np.abs((data - data[-1]) / data[-1])

            # Evaluate what is the lower converged data
            i = np.argwhere(error_to_last <= p['conv_thr'])
            index_converged = index[i]
            

            if p['parallelization'] == 'low':
                if len(index_converged > 1):
                    self.store_to_db(index_converged[0])
                    stop_convergence = False
            
            elif p['parallelization'] == 'high':
                self.store_to_db(index_converged[0])
                stop_convergence = True

            return stop_convergence

        def store_to_db(self, index):

            nav = Navigator(db_file=p['db_file'], high_level=p['high_level'])
            out_dict = nav.find_data(collection=p['functional'] + '.slab_data', 
                                     filter={'mpid': p['mp_id'], 'miller': p['miller']})
            
            # Extract the data to be saved elsewhere
            entry = ['thickness', 'data_' + str(index), 'calc_output']
            store_dict = get_one_info_from_dict(out_dict, entry)

            nav.update_data(collection=p['functional'] + '.slab_data', 
                            filter={'mpid': p['mp_id'], 'miller': p['miller']},
                            new_values={'$set': {'calc_output': store_dict}})

        def call_recursion(self, fw_spec, p):
            
            from triboflow.workflows.slabs_wfs import SlabWF

            # Select the correct function to call the workflow 
            if p['conv_kind'] == 'surfene':
                generate_wf = SlabWF.conv_slabthick_surfene
            else:
                pass
            
            # Generate the workflow for the detour
            wf = generate_wf(structure=p['structure'], mp_id=p['mp_id'], 
                             miller=p['miller'], functional=p['functional'], 
                             comp_params=p['comp_params'], spec=fw_spec, 
                             db_file=p['db_file'], low_level=p['low_level'], 
                             high_level=p['high_level'], relax_type=p['relax_type'], 
                             thick_min=p['thick_min'], thick_max=p['thick_max'], 
                             thick_incr=p['thick_incr'], vacuum=p['vacuum'], 
                             in_unit_planes=p['in_unit_planes'], ext_index=p['ext_index'], 
                             conv_thr=p['conv_thr'], parallelization=p['parallelization'], 
                             recursion=True, cluster_params=p['cluster_params'])
            
            return wf
