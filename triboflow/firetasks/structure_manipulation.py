""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""
import monty
import numpy as np
from uuid import uuid4
from pprint import pprint, pformat

from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from fireworks import FWAction, FiretaskBase, Firework, Workflow, FileWriteTask
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk

from triboflow.workflows.base import dynamic_relax_swf
from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.vasp_tools import get_custom_vasp_relax_settings
from triboflow.utils.structure_manipulation import (
    interface_name, slab_from_structure, recenter_aligned_slabs, 
    stack_aligned_slabs, transfer_average_magmoms, clean_up_site_properties)
from triboflow.utils.file_manipulation import copy_output_files
from triboflow.phys.interface_matcher import MatchInterface




@explicit_serialize
class FT_StartPreRelax(FiretaskBase):
    """ Start a subworkflow as a detour to relax the cell shape and positions
    of a primitive structure depending on lattice parameters, and then move
    the optimized primitive structure to the high level database.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    encut : float, optional
        Energy cutoff for the relaxation run. Defaults to 1000.
    k_dens : int, optional
        kpoint density in 1/Angstrom. Defaults to a (quite high) 20.
    high_level_db : str, optional
        Name of the high level database the structure should be queried in
        and later the results written to. Defaults to 'triboflow'.
    Returns
    -------
        Starts a subworkflow as a detour to the current workflow.
    """
    _fw_name = 'Start a cell shape relaxation'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file', 'encut', 'k_dens', 'high_level_db']

    def run_task(self, fw_spec):
        mp_id = self.get('mp_id')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        # Querying the structure from the high level database.
        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        data = nav_structure.get_bulk_from_db(
            mp_id=mp_id,
            functional=functional)

        prim_struct = data.get('primitive_structure')
        if not prim_struct:
            struct = data.get('structure_fromMP')
            if not struct:
                raise LookupError('No structure found in the database that can '
                                  'be used as input for cell shape relaxation.')
            struct = Structure.from_dict(struct)
            prim_struct = SpacegroupAnalyzer(struct).get_primitive_standard_structure()
            prim_struct = transfer_average_magmoms(struct, prim_struct)
        else:
            prim_struct = Structure.from_dict(prim_struct)

        a = np.round(prim_struct.lattice.a, 6)
        b = np.round(prim_struct.lattice.b, 6)
        c = np.round(prim_struct.lattice.c, 6)

        if data.get('pre_relaxed') or ((a == c) and (b == c)):
            return FWAction(update_spec=fw_spec)
        else:

            # Querying the computational parameters from the high level database
            # and updating with the optional inputs
            comp_params = data.get("comp_parameters")
            encut = self.get('encut', 1000)
            k_dens = self.get('k_dens', 15)
            comp_params.update({"encut": encut, "k_dens": k_dens})

            tag = "CellShapeRelax-{}".format(str(uuid4()))
            vis = get_custom_vasp_relax_settings(prim_struct, comp_params,
                                                 'bulk_pos_shape_relax')
            RelaxWF = dynamic_relax_swf([[prim_struct, vis, tag]])

            MoveResultsFW = Firework([FT_UpdatePrimStruct(functional=functional,
                                                          tag=tag, flag=mp_id,
                                                          high_level_db=hl_db)],
                                     name='Move pre-relaxed structure for {}'.
                                     format(prim_struct.formula))
            MoveResultsWF = Workflow([MoveResultsFW])

            RelaxWF.append_wf(MoveResultsWF, RelaxWF.leaf_fw_ids)

            return FWAction(detours=RelaxWF, update_spec=fw_spec)


@explicit_serialize
class FT_UpdatePrimStruct(FiretaskBase):
    """ Update the primitive structure in the high level database.

    Parameters
    ----------
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    flag : str
        An identifyer to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise use something
        unique which you can find again.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    Returns
    -------
        Moves the optimized primitive structure from the low level
        database to the high level database.
    """
    _fw_name = 'Update primitive structure in the high level DB'
    required_params = ['functional', 'tag', 'flag']
    optional_params = ['db_file', 'high_level_db']

    def run_task(self, fw_spec):
        functional = self.get('functional')
        tag = self.get('tag')
        flag = self.get('flag')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        nav = Navigator(db_file=db_file)
        calc = nav.find_data('tasks', {'task_label': tag})
        out = calc['output']

        struct_dict = {'primitive_structure': out['structure'],
                       'pre_relaxed': True}

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        nav_high.update_data(
                    collection=functional+'.bulk_data',
                    fltr={'mpid': flag},
                    new_values={'$set': struct_dict},
                    upsert=False)

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartSlabRelaxSWF(FiretaskBase):
    """Start a subworkflow as a detour to make a slab and relax it.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    slab_struct_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'unrelaxed_slab'.
    relax_type : str
        Relaxation type for the get_custom_vasp_relax_settings helper_function.
    bulk_struct_name : str, optional
        Name of the bulk structure in the bulk database (material is
        identified by mp_id, but there might be different structures of the
        same material.) Defaults to 'structure_equiVol'.
    slab_out_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'relaxed_slab'.
        
    Returns
    -------
        Starts a new subworkflow as a detour to the current workflow.
    """
    
    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'slab_struct_name', 'relax_type',
                       'bulk_struct_name', 'slab_out_name', 'high_level_db']

    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import make_and_relax_slab_swf

        mp_id = self.get('mp_id')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
        else:
            miller = self['miller']

        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        relax_type = self.get('relax_type', 'slab_pos_relax')
        bulk_name = self.get('bulk_struct_name', 'structure_equiVol')
        slab_out_name = self.get('slab_out_name', 'relaxed_slab')
        hl_db = self.get('high_level_db', True)
        
        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)

        data = nav_structure.get_bulk_from_db(
            mp_id=mp_id,
            functional=functional)
        bulk_struct = Structure.from_dict(data[bulk_name])
        
        slab_data = nav_structure.get_slab_from_db(
            mp_id=mp_id,
            functional=functional,
            miller=miller)
        comp_params = slab_data.get('comp_parameters')
        min_thickness = slab_data.get('min_thickness', 10)
        min_vacuum = slab_data.get('min_vacuum', 25)
        
        WF = make_and_relax_slab_swf(bulk_structure=bulk_struct,
                                     miller_index=miller,
                                     flag=mp_id,
                                     comp_parameters=comp_params,
                                     functional=functional,
                                     min_thickness=min_thickness,
                                     min_vacuum=min_vacuum,
                                     relax_type=relax_type,
                                     slab_struct_name=slab_name,
                                     out_struct_name=slab_out_name,
                                     print_help=False)
        
        return FWAction(detours=WF)
    
@explicit_serialize
class FT_GetRelaxedSlab(FiretaskBase):
    """Get the relaxed structure from the DB, and put a Slab into the high-level DB.

    Parameters
    ----------
    flag : str
        An identifyer to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise use something
        unique which you can find again.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    struct_out_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'relaxed_slab'.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The user name on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.
    
    Returns
    -------
        Pymatgen Slab into the high-level DB.
    """
    
    required_params = ['flag', 'miller', 'functional', 'tag']
    optional_params = ['db_file', 'struct_out_name', 'file_output', 'high_level_db',
                       'output_dir', 'remote_copy', 'server', 'user', 'port']

    def run_task(self, fw_spec):

        flag = self.get('flag')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
        else:
            miller = self['miller']

        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        hl_db = self.get('high_level_db', True)
        out_name = self.get('struct_out_name', 'relaxed_slab')
        file_output = self.get('file_output', False)
        output_dir = self.get('output_dir', None)
        remote_copy = self.get('remote_copy', False)
        server = self.get('server', None)
        user = self.get('user', None)
        port = self.get('port', None)
        
        # Check if a relaxed slab is already in the DB entry
        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        slab_data = nav_structure.get_slab_from_db(
            mp_id=flag,
            functional=functional,
            miller=miller)
        
        if out_name not in slab_data:
            # Get results from OptimizeFW
            nav = Navigator(db_file=db_file)
            vasp_calc = nav.find_data(
                collection='tasks', 
                fltr={'task_label': self['tag']})
            relaxed_slab = Structure.from_dict(vasp_calc['output']['structure'])
            slab = Slab(relaxed_slab.lattice,
                        relaxed_slab.species_and_occu,
                        relaxed_slab.frac_coords,
                        miller,
                        Structure.from_sites(relaxed_slab, to_unit_cell=True),
                        shift=0,
                        scale_factor=[[1,0,0], [0,1,0], [0,0,1]],
                        site_properties=relaxed_slab.site_properties)
            
            nav_high = Navigator(db_file=db_file, high_level=hl_db)
            nav_high.update_data(
                collection=functional+'.slab_data',
                fltr={'mpid': flag, 'miller': miller},
                new_values={'$set': {out_name: slab.as_dict()}})
        else:
            nav = Navigator(db_file=db_file)
            vasp_calc = nav.find_data(
                collection='tasks',
                fltr={'task_label': self['tag']})

            if  vasp_calc:
                relaxed_slab = Structure.from_dict(vasp_calc['output']['structure'])
                slab = Slab(relaxed_slab.lattice,
                        relaxed_slab.species_and_occu,
                        relaxed_slab.frac_coords,
                        miller,
                        Structure.from_sites(relaxed_slab, to_unit_cell=True),
                        shift=0,
                        scale_factor=[[1,0,0], [0,1,0], [0,0,1]],
                        site_properties=relaxed_slab.site_properties)
                print('')
                print(' A slab with the selected output name already exists in the DB.'
                      ' It will not be overwritten with the new relaxed slab.\n'
                      ' If needed you can update the data manually.')
                print('')
            else:
                slab = Slab.from_dict(slab_data[out_name])
                print('')
                print(' A slab with the selected output name already exists in the DB.'
                      ' No new slab has been relaxed.\n')

        # screen output:
        print('')
        print('Relaxed output structure as pymatgen.surface.Slab dictionary:')
        pprint(slab.as_dict())
        print('')
        
        # handle file output:
        if file_output:
            poscar_str = Poscar(slab).get_string()
            poscar_name = flag+'_Relaxed_slab_POSCAR.vasp'
            slab_name = flag+'_Relaxed_slab_dict.txt'
            write_FT = FileWriteTask(
                files_to_write=[{'filename': poscar_name,
                                 'contents': poscar_str},
                                {'filename': slab_name,
                                 'contents': pformat(slab.as_dict())}])
            copy_FT = copy_output_files(file_list=[poscar_name, slab_name],
                                        output_dir=output_dir,
                                        remote_copy=remote_copy,
                                        server=server,
                                        user=user,
                                        port=port)
            FW = Firework([write_FT, copy_FT],
                          name='Copy SlabRelax SWF results')
            WF = Workflow.from_Firework(FW, name='Copy SlabRelax SWF results')

            return FWAction(update_spec=fw_spec, detours=WF)
        else:  
            return FWAction(update_spec=fw_spec)

            
@explicit_serialize
class FT_StartSlabRelax(FiretaskBase):
    """Fetch a Slab from the DB and relax it using an atomate OptimizeFW.

    Parameters
    ----------
    flag : str
        An identifyer to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise use something
        unique which you can find again.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.
    slab_struct_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'unrelaxed_slab'.
    relax_type : str
        Relaxation type for the get_custom_vasp_relax_settings helper_function.
    
    Returns
    -------
        Pymatgen Slab into the high-level DB.
    """
    
    required_params = ['flag', 'miller', 'functional', 'tag']
    optional_params = ['db_file', 'comp_parameters', 'slab_struct_name',
                       'relax_type', 'high_level_db']

    def run_task(self, fw_spec):

        flag = self.get('flag')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
            miller_str = self['miller']
        else:
            miller = self['miller']
            miller_str = ''.join(str(s) for s in self['miller'])

        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        comp_params = self.get('comp_parameters', {})
        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        tag = self.get('tag')
        relax_type = self.get('relax_type', 'slab_pos_relax')
        hl_db = self.get('high_level_db', True)
        
        nav_structure = StructureNavigator(
            db_file=db_file, 
            high_level=hl_db)
        slab_data = nav_structure.get_slab_from_db(
            mp_id=flag,
            functional=functional,
            miller=miller)
        
        slab_to_relax = Slab.from_dict(slab_data.get(slab_name))
        formula = slab_to_relax.composition.reduced_formula
        
        # Check if a relaxed slab is already in the DB entry
        if 'relaxed_slab' not in slab_data:
            vis = get_custom_vasp_relax_settings(slab_to_relax, comp_params,
                                                 relax_type)
            inputs = [[slab_to_relax, vis, tag]]
            wf_name = formula+miller_str+'_'+relax_type
            WF = dynamic_relax_swf(inputs_list=inputs,
                                   wf_name=wf_name)
                        
            return FWAction(detours=WF)

@explicit_serialize
class FT_MakeSlabInDB(FiretaskBase):
    """Makes a slab with certain orientation out of a bulk structure.

    The slab has a defined thickness and orientation.

    Parameters
    ----------
    bulk_structure : pymatgen.core.structure.Structure
        Bulk structure that is used to construct the slab out of.
    miller : list of int or str
        Miller indices of the slab to make.
    flag : str
        An identifyer to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise use something
        unique which you can find again.
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    slab_struct_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'unrelaxed_slab'.
    min_thickness : float, optional
        Minimal thickness of the unrelaxed slab in Angstrom. The default is 10.0.
    min_vacuum : float, optional
        Minimal thickness of the vacuum layer in Angstrom. The default is 25.0.
    
    Returns
    -------
        Pymatgen Slab into the high-level DB.
    """
    
    _fw_name = 'Make slab from bulk structure'
    required_params = ['bulk_structure', 'miller', 'flag', 'functional']
    optional_params = ['db_file', 'slab_struct_name', 'min_thickness',
                       'min_vacuum', 'high_level_db']
    
    def run_task(self, fw_spec):
        
        bulk_prim = self.get('bulk_structure')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
        else:
            miller = self['miller']
        flag = self.get('flag')
        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        min_thickness = self.get('min_thickness', 10)
        min_vacuum = self.get('min_vacuum', 25)
        hl_db = self.get('high_level_db', True)
        
        # nav_structure = StructureNavigator(
        #     db_file=db_file,
        #     high_level=hl_db)
        # slab_data = nav_structure.get_slab_from_db(
        #     mp_id=flag, 
        #     functional=functional,
        #     miller=miller)
        
        bulk_conv = SpacegroupAnalyzer(bulk_prim).get_conventional_standard_structure()
        bulk_conv = transfer_average_magmoms(bulk_prim, bulk_conv)
        
        SG = SlabGenerator(initial_structure=bulk_conv,
                           miller_index=miller,
                           center_slab=True,
                           primitive=True,
                           lll_reduce=True,
                           max_normal_search=max([abs(l) for l in miller]),
                           min_slab_size=min_thickness,
                           min_vacuum_size=min_vacuum)
        
# =============================================================================
# Here we need to do more Work! Now the first slab in the list is taken,
# which should be fine for monoatomic slabs, or if the methods only finds
# one slab which has no broken bonds, but in general, several slabs should be
# investigated for the lowest surface energy, and this one should be taken.
# =============================================================================
        slab = SG.get_slabs(bonds=None, ftol=0.1, tol=0.1, max_broken_bonds=0,
                            symmetrize=False, repair=False)[0]
        
        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        slab_dict = monty.json.jsanitize(slab.as_dict(), allow_bson=True)
        nav_high.update_data(
            collection=functional+'.slab_data',
            fltr={'mpid': flag, 'miller': miller},
            new_values={'$set': {slab_name: slab_dict}},
            upsert=True)

        return


@explicit_serialize
class FT_MakeHeteroStructure(FiretaskBase):
    """Matches two slab systems to form a heterogeneous interface.

    If the match fails, the accuracy criteria are relaxed in steps of 5% until
    a match is found.
    
    Parameters
    ----------
    mp_id_1 : str
        ID number from the MP to identify material 1.
    miller_1 : list of int or str
        Miller indices for identification of the slab of material 1. Either
        single str e.g. '111', or a list of int e.g. [1, 1, 1].
    mp_id_2 : str
        ID number from the MP for material 2.
    miller_2 : list of int or str
        Miller indices for identification of the slab of material 2. Either
        single str e.g. '111', or a list of int e.g. [1, 1, 1].
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.

    Returns
    -------
    Matched interface structure and bottom and top slabs in interface_data DB.
    """
    
    _fw_name = 'Make Hetero Structure'
    required_params = ['mp_id_1', 'miller_1', 'mp_id_2', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'high_level_db']
    
    def run_task(self, fw_spec):

        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        if type(self['miller_1']) == str:
            miller_1 = [int(k) for k in list(self['miller_1'])]
        else:
            miller_1 = self['miller_1']
        if type(self['miller_2']) == str:
            miller_2 = [int(k) for k in list(self['miller_2'])]
        else:
            miller_2 = self['miller_2']
        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
            
        hl_db = self.get('high_level_db', True)

        nav_high = Navigator(db_file=db_file, high_level=hl_db)
        inter_name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)
        inter_data = nav_high.find_data(
            collection=functional+'.interface_data',
            fltr={'name': inter_name})
        
        inter_params = inter_data['interface_parameters']

        if not inter_data.get('unrelaxed_structure'):
            
            nav_structure = StructureNavigator(
                db_file=db_file, 
                high_level=hl_db)
            slab_1_dict = nav_structure.get_slab_from_db(
                mp_id=mp_id_1,
                functional=functional,
                miller=miller_1)
            slab_2_dict = nav_structure.get_slab_from_db(
                mp_id=mp_id_2,
                functional=functional,
                miller=miller_2)

            slab_1 = Slab.from_dict(slab_1_dict['relaxed_slab'])
            slab_2 = Slab.from_dict(slab_2_dict['relaxed_slab'])
            
            bulk_dat_1 = nav_structure.get_bulk_from_db(mp_id_1, functional)
            bulk_dat_2 = nav_structure.get_bulk_from_db(mp_id_2, functional)
            
            bm_1 = bulk_dat_1['bulk_moduls']
            bm_2 = bulk_dat_2['bulk_moduls']
# =============================================================================
# Running into crashes for max_angle_diff > ~1.5 for MPInterfaces 2020.6.19,
# at least for certain interfaces. A match is found, but than there is a 
# LinAlgError("Singular matrix") error in forming the matched slabs?
# The following lines ensures that max_angle_diff > 1.5. This is not a great
# solution obviously. I also changed the default in FT_CheckInterfaceParamDict
#
# Ran into this problem again only more severe in November 2020. For Ag111 and
# Pt111 matching the error appears for max_angle_diff > ~0.4!
# Please see issue #25 on gitlab.
# =============================================================================
            # if inter_params['max_angle_diff'] > 0.4:
            #     inter_params['max_angle_diff'] = 0.4
        
            MI = MatchInterface(slab_1=slab_1,
                                slab_2=slab_2,
                                strain_weight_1=bm_1,
                                strain_weight_2=bm_2,
                                **inter_params)
            top_aligned, bottom_aligned = MI.get_centered_slabs()
        
            if top_aligned and bottom_aligned:
                
                interface = MI.get_interface()
                
                inter_dict = interface.as_dict()
                bottom_dict = bottom_aligned.as_dict()
                top_dict = top_aligned.as_dict()

                nav_high.update_data(
                    collection=functional+'.interface_data',
                    fltr={'name': inter_name},
                    new_values={'$set': {'unrelaxed_structure': inter_dict,
                                         'bottom_aligned': bottom_dict,
                                         'top_aligned': top_dict}})
            
            else:

                f = 1.05 # factor with which to multiply the missmatch criteria
            
                new_params = {
                    'max_area': inter_params['max_area'] * f,
                    'max_mismatch': inter_params['max_mismatch'] * f,
                    'max_angle_diff': inter_params['max_angle_diff'] * f,
                    'r1r2_tol': inter_params['r1r2_tol'] * f,
                    'interface_distance': inter_params['interface_distance']}
            
                if not inter_data.get('original_interface_params'):

                    nav_high.update_data(
                        collection=functional+'.interface_data',
                        fltr={'name': interface_name},
                        new_values={'$set': 
                                       {'original_interface_params': 
                                            inter_params}})

                nav_high.update_data(
                    collection=functional+'.interface_data',
                    fltr={'name': interface_name},
                    new_values={'$set': {'interface_parameters': new_params}})

                new_fw = Firework(FT_MakeHeteroStructure(mp_id_1=mp_id_1,
                                                         miller_1=miller_1,
                                                         mp_id_2=mp_id_2,
                                                         miller_2=miller_2,
                                                         db_file=db_file))
                
                return FWAction(detours=new_fw)



@explicit_serialize
class FT_AddSelectiveDynamics(FiretaskBase):
    """Write a POSCAR with selective dynamics added.
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure from which the POSCAR file is written
    selective_dynamics_array: list, optional
        Optional Nx3 array of booleans given as a list of lists. N is equal to
        the number of sites in the structure.
    ----------
    
    Returns
    -------
    POSCAR: file
        A POSCAR file is written in the current directory.
    
    """
    
    _fw_name = 'Add selective dynamics'
    required_params = ['structure']
    optional_params = ['selective_dynamics_array']

    def run_task(self, fw_spec):
        struct = self['structure']
        
        if 'selective_dynamics_array' in self:
            sd_array = self['selective_dynamics_array']
            if len(sd_array) != len (struct.sites):
                raise SystemExit('Length of provided selective_dynamics array'
                                 ' is not equal to the number of sites in the'
                                 ' structure!')
        else:
            sd_array = []
            for i in range(len(struct.sites)):
                sd_array.append([False, False, True])
            #sd_array = np.asarray(sd_array)
        
        poscar = Poscar(struct, selective_dynamics=sd_array)
        poscar.write_file('POSCAR')
        
        spec = fw_spec

        return FWAction(update_spec = spec)