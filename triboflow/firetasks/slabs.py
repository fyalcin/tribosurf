#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:09:03 2021

Firetasks to converge the slab thickness either calculating the surface 
energy or the lattice parameter for structures with increasing number of 
atomic layers.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

import os
import monty
from uuid import uuid4

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction
from fireworks import Workflow
from atomate.utils.utils import env_chk
from atomate.vasp.powerups import add_modify_incar
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator

from triboflow.utils.database import StructureNavigator
from triboflow.utils.errors import SlabOptThickError, GenerateSlabsError, RelaxStructureError, ReadParamsError
from triboflow.tasks.io import read_json
from triboflow.utils.vasp_tools import GetCustomVaspRelaxSettings
from atomate.vasp.fireworks.core import OptimizeFW, ScanOptimizeFW

currentdir = os.path.dirname(__file__)

# ============================================================================
# Firetasks
# ============================================================================

@explicit_serialize
class FT_SlabOptThick(FiretaskBase):
    """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

    Start a subworkflow as a detour to calculate the optimal thickness generated
    from a provided bulk structure and with a certain orientation.
    The thickness can be converged either by evaluating how the surface energy
    or the lattice parameters changes with the number of atomic planes.
    At the moment only the surface energy convergence is implemented.

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
        surface energy will be saved too if convo_kind is 'surfene'.
        The default is 'triboflow'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.

    convo_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene', 
        'alat'. The latter is not implemented yet. The default is 'surfene'.
    
    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point.
        The default is 4.
        
    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If the convergence
        is not reached this value is assumed to be the optimal one.
        The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step.
        The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.

    bulk_name : str or list None, optional
        Name of the custom bulk dictionary, to be retrieved from the high level
        database and to be used to build the slabs. Bulks are identified by 
        mp_id and functional but there might be different structures of the
        same material. The default is "structure_fromMP".

    slab_name : str or None, optional
        Where to search for the information about the optimal thickness in the
        slab dictionary in the high level database. The default is None.

    """
    
    _fw_name = 'Start a subworkflow to converge slab thickness'

    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'low_level', 'high_level', 'convo_kind',
                       'relax_type', 'thick_min', 'thick_max', 'thick_incr',
                       'vacuum', 'in_unit_planes', 'bulk_name', 'slab_name']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file=dfl, default_key="SlabOptThick")
        
        # Retrieve the bulk information from the high level DB
        nav_struct = StructureNavigator(p['db_file'], p['high_level'])
        slab = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                           functional=p['functional'], 
                                           miller=p['miller'])
        
        # If data is saved elsewhere from standard position
        if p['slab_name'] is not None:
            slab = one_info_from_struct_dict(slab, p['slab_name'])

        # Start a subworkflow to converge the thickness if not already done
        stop_convergence = slab.get('opt_thickness', None)

        if not stop_convergence:
            
            # Retrieve the desired bulk structure
            bulk = nav_struct.get_slab_from_db(mp_id=p['mp_id'], 
                                               functional=p['functional'])
            bulk = one_info_from_struct_dict(bulk, p['bulk_name'])

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

        from triboflow.firetasks.slabs_wfs import SlabWF

        if p['convo_kind'] == 'surfene':
            generate_wf = SlabWF.conv_slabthick_surfene
        elif p['convo_kind'] == 'alat':
            generate_wf = SlabWF.conv_slabthick_alat
        else:
            raise SlabOptThickError("Wrong input argument for convo_kind. "
                                    "Allowed options: 'surfene', 'alat'.")

        wf = generate_wf(structure=structure, mp_id=p['mp_id'], 
                         miller=p['miller'], functional=p['functional'], 
                         comp_params=p['comp_params'], db_file=p['db_file'],
                         low_level=p['low_level'], high_level=p['high_level'],
                         relax_type=p['relax_type'], thick_min=p['thick_min'], 
                         thick_max=p['thick_max'], thick_incr=p['thick_incr'], 
                         vacuum=p['vacuum'], in_unit_planes=p['in_unit_planes'])

        return wf

@explicit_serialize
class FT_StartThickConvo(FiretaskBase):
    """ Author: Gabriele Losi; Copyright 2021, Prof. M.C. Righi, UniBO.

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
        
    database : str or None, optional
        Name of the table where the data will be stored, within db_file. If 
        nothing is passed the Firework database is used. The default is None.

    convo_kind : str, optional
        Type of convergence to be performed. Allowed values are: 'surfene', 
        'alat'. The latter is not implemented yet. The default is 'surfene'.

    relax_type : str, optional
        The type of relaxation to be performed during the simulation, to be feed
        to GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.

    comp_params : dict, optional
        Computational parameters for the VASP simulations. If not set, default 
        parameters will be used instead. The default is {}.
    
    thick_min : int, optional
        Number of atomic layers for the slab to be used as starting point.
        The default is 4.
        
    thick_max : int, optional
        Maximum number of allowed atomic layers for the slab. If the convergence
        is not reached this value is assumed to be the optimal one.
        The default is 12.

    thick_incr : int, optional
        The incremental number of atomic layers to be added at each step.
        The default is 2.

    vacuum : int or float, optional
        Vacuum to be used for creating in the slabs cells. The default is 10.

    in_unit_planes : bool, optional
        Decide if thick_min, thick_max, thick_incr, and vacuum are expressed in
        units of number of atomic layers or Angstrom. The default is True.

    ext_index : int, optional
        Use the ext_index element from SlabGenerator.get_slabs as a slab.
        The default is 0.

    slab_name : str or None, optional
        Custom name to store the slab dictionary in the high level database at
        the end of the convergence procedure. The default is None.

    cluster_params : dict, optional
        Dictionary containing cluster-related options to run efficiently the
        VASP simulations on a cluster. The default is {}.

    """
    
    _fw_name = 'Start the slab thickness convergence'
    required_params = ['structure', 'mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'database', 'convo_kind', 'relax_type',
                       'comp_params', 'thick_min', 'thick_max', 'thick_incr',
                       'vacuum', 'in_unit_planes', 'ext_index', 'slab_name', 
                       'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """
        from triboflow.firetasks.surfene_wfs import SurfEneWF

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file=dfl, default_key="StartThickConvo")

        if p['convo_kind'] == 'surfene':
            wf = SurfEneWF.surface_energy(structure=p['structure'],
                                          mp_id=p['mp_id'], 
                                          miller=p['miller'], 
                                          functional=p['functional'],
                                          comp_params=p['comp_params'],
                                          db_file=p['db_file'], 
                                          database=p['database'],
                                          thick_start=p['thick_start'], 
                                          thick_incr=p['thick_incr'], 
                                          nsteps=p['nsteps'],
                                          vacuum=p['vacuum'],
                                          relax_type=p['relax_type'],
                                          ext_index=p['ext_index'], 
                                          slab_name=p['slab_name'],
                                          cluster_params=p['cluster_params'])
        else:
            raise SystemExit('Lattice parameter convo not yet implemented')
        
        return FWAction(detours=wf, update_spec=fw_spec)

@explicit_serialize
class FT_EndThickConvo(FiretaskBase):
    """
    Firetask description...
    
    """

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 
    pass

@explicit_serialize
class FT_GenerateSlabs(FiretaskBase):
    """
     Generate a slab or a list of slabs out of a given structure. Parameters 
     hat are taken into account to generate the possible different slabs are: 
     miller, thickness, vacuum, slab_name. The slabs are generated with 
     SlabGenerator and stored in the database.
    
    """

    required_params = ['structure', 'mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'database', 'thickness', 'vacuum', 
                       'symmetrize', 'ext_index', 'tag', 'slab_name']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file=dfl, default_key="GenerateSlabs")

        GenerateSlabsError.check_slabname_thickness(p['thickness'],
                                                    p['slab_name'])

        # Generate the slabs for each structure passed as input
        slabs = []
        thickness = list(p['thickness'])
        slab_name = list(p['slab_name'])

        # # Create the bulk
        # slabgen = SlabGenerator(initial_structure = p['structure'],
        #                         miller_index = p['miller'],
        #                         primitive=False,
        #                         lll_reduce=True,
        #                         in_unit_planes=True,  # Fundamental
        #                         min_slab_size=thickness[0],
        #                         min_vacuum_size=0)

        for thk in thickness:

            if thk == 0:
                pass
                #bulk = slabgen.oriented_unit_cell()
                
            slabgen = SlabGenerator(initial_structure = p['structure'],
                                    miller_index = p['miller'],
                                    center_slab=True,
                                    primitive=False,
                                    lll_reduce=True,
                                    in_unit_planes=True,  # Fundamental
                                    min_slab_size=thk,
                                    min_vacuum_size=p['vacuum'])
            
            s = slabgen.get_slabs(bonds=None, ftol=0.1, tol=0.1, 
                                  max_broken_bonds=0, symmetrize=p['symmetrize'],
                                  repair=False)
            # If exit index is True (default), return the first element.
            # TODO: Generalize this approach by evaluating the different slabs
            # that are returned and select the optimal one.
            if p['ext_index']:
                s = s[0]

            slabs.append(s)

        # Add the slabs to the "database" DB within the db_file path
        nav_struct = StructureNavigator(p['db_file'], p['database'])
        
        # Store unrelaxed data in the Database
        for slab, name in zip(slabs, slab_name):
            slab_dict = monty.json.jsanitize(s.as_dict(), allow_bson=True)
            nav_struct.add_slab_to_db(structure=slab, 
                                      mp_id=p['mp_id'],
                                      functional=p['functional'],
                                      miller=p['miller'],
                                      struct_name=name)
        
        # coll.update_one({'mpid': flag, 'miller': miller},
        #                 {'$set': {'unrelaxed_slab': slab_dict}},
        #                 upsert=True)
        
        return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_RelaxStructure(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['mp_id', 'functional', 'struct_kind']
    optional_params = ['comp_params', 'miller', 'name', 'db_file', 'database',
                       'relax_type']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file=dfl, default_key="RelaxStructure")

        # Get the structure from the Database with a calculation tag
        structure, tag = self.extract_data_from_db()
        if p['name'] is not None:
            structure = structure['name']

        # Set the simulation if the calculation is not already at convergence
        if 'relaxed' or 'relaxed_slab' not in structure:
            fw = self.set_calculation(structure, tag, dfl)

            # Set the Relaxation Workflow for making a detour.
            # Use add_modify_incar powerup to add KPAR and NCORE settings based on 
            # env_chk in my_fworker.yaml.
            wf_name = tag.split('_')[0] + p['relax_type']
            opt_wf = add_modify_incar(Workflow([fw], name=wf_name))

            return FWAction(detours=opt_wf)
    
    def extract_data_from_db(self):

        RelaxStructureError.check_struct_kind(self.p['struct_kind'])

        # Call the navigator for retrieving structure
        nav_struct = StructureNavigator(db_file=self.p['db_file'],
                                        high_level=self.p['database'])

        # Extract data from database
        if self.p['struct_kind'] == 'bulk':
            structure = nav_struct.get_bulk_from_db(self.p['mp_id'], 
                                                    self.p['functional'],
                                                    warning=True)
        elif self.p['struct_kind'] == 'slab':
            structure = nav_struct.get_slab_from_db(self.p['mp_id'], 
                                                    self.p['functional'], 
                                                    self.p['miller'], 
                                                    warning=True)
        elif self.p['struct_kind'] == 'interface':
            structure = nav_struct.get_interface_from_db(self.p['name'],
                                                         self.p['functional'],
                                                         warning=True)
        
        RelaxStructureError.is_data(structure=structure, mp_id=self.p['mp_id'], 
                                    functional=self.p['functional'], 
                                    struct_kind=self.p['struct_kind'])
        
        # Get a tag for the calculation
        formula = structure.composition.reduced_formula
        if self.p['miller'] is not None:
            miller_str = ''.join(str(s) for s in self.p['miller'])
            tag = formula + miller_str + '_' + str(uuid4())
        else:
            tag = formula + '_' + str(uuid4())

        return structure, tag

    def set_calculation(self, structure, tag, dfl):

        # Check the computational parameters
        comp_params = self.p['comp_params']
        if not comp_params:
            defaults = read_json(dfl)
            comp_params = defaults['comp_params']

        # Set options for vasp
        vis = GetCustomVaspRelaxSettings(structure, comp_params,
                                         self.p['relax_type'])
        
        # Create the Firework to run the simulation
        if self.p['functional'] == 'SCAN':
            fw = ScanOptimizeFW(structure=structure, name=tag, vasp_input_set=vis)
        else:
            fw = OptimizeFW(structure, name=tag, vasp_input_set=vis,
                            half_kpts_first_relax=True)

        return fw

@explicit_serialize
class FT_PutStructInDB(FiretaskBase):
    """
    Firetask description...
    
    """

    required_params = ['mp_id', 'functional', 'tag', 'struct_kind']
    optional_params = ['miller', 'name', 'db_file', 'cluster_params']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, required_params, optional_params,
                                default_file=dfl, default_key="PutStructInDB")

        # Set missing cluster params
        self.check_cluster_params(self, dfl)

        # Get the structure from the Database with a calculation tag
        structure, tag = self.extract_data_from_db()
        if p['name'] is not None:
            structure = structure['name']

    def check_cluster_params(self, dfl):
        
        defaults = read_json(dfl)
        cluster_params = self.p['cluster_params']

        # Set default parameters if key is absent
        if not cluster_params:
            cluster_params = defaults['cluster_params']
        
        else:
            for key, value in defaults.items():
                if not key in cluster_params.keys():
                    cluster_params[key] = value
        
        self.p['cluster_params'] = cluster_params

# ============================================================================
# Functions
# ============================================================================

def read_runtask_params(obj, fw_spec, required_params, optional_params,
                        default_file, default_key):
    """
    General function to read the required and optional parameters of a 
    firetask instance. If 

    Parameters
    ----------
    obj : [type]
        [description]
    fw_spec : [type]
        [description]
    required_params : [type]
        [description]
    optional_params : [type]
        [description]
    default_file : [type]
        [description]
    default_key : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    defaults = read_json(default_file)
    params = {}

    # Read required and optional parameters
    for key in required_params:
        params[key] = obj[key]
    for key in optional_params:
        params[key] = obj.get(key, defaults[default_key][key])

    # Clean db_file
    if 'db_file' in params.keys():
        if params['db_file'] is None:
            params['db_file'] = env_chk('>>db_file<<', fw_spec)

    # Clean miller index
    if 'miller' in params.keys():
        if isinstance(params['miller'], str):
            miller = [int(k) for k in list(params['miller'])]

    return params

def one_info_from_struct_dict(struct_dict, name):

    # Simply read a dictionary key
    if isinstance(name, str):
        info = struct_dict[name]
    
    # You can have multiple innested keys
    elif isinstance(name, list):
        if all([isinstance(x, str) for x in name]):

            info = struct_dict.copy()
            for n in name:
                info = info[n]  # Read one key after the other
    
    else:
        ReadParamsError('Error in reading struct_dict, name is wrong.')

    return info

def multiple_info_from_struct_dict(struct_dict, name):

    # Extract many info at the same time 
    if isinstance(name, list):
        if all([isinstance(x, list) for x in name]):
            info = []
            for n in name:
                info.append(one_info_from_struct_dict(struct_dict, n))

        else:
            info = one_info_from_struct_dict(struct_dict, name)
    
    return info


# Small test
if __name__ == '__main__':
    FT_SlabOptThick(mp_id='mp-100', functional='PBE', miller=[1,0,0])
