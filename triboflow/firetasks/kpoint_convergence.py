"""Firetasks for converging the kpoints density using absolute energy differences
 in the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from datetime import datetime
from pprint import pformat

from pymatgen.core.structure import Structure
from fireworks import FWAction, FiretaskBase, Firework, Workflow, FileWriteTask
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.powerups import add_modify_incar

from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.check_convergence import is_list_converged
from triboflow.utils.vasp_tools import (get_custom_vasp_static_settings,
    MeshFromDensity)
from triboflow.utils.file_manipulation import copy_output_files


@explicit_serialize
class FT_StartKPointConvo(FiretaskBase):
    """ Starts a kpoint density convergence subworkflow.

    THIS IS A DEPRECIATED FIRETASK THAT HAS BEEN REPLACED BY:
        FT_StartConvo in triboflow.firetasks.convergence
    Its function is to perform a kpoint density convergence of a given
    material identified by its MPID with the given functional.

    Parameters
    ----------
    mp_id : str
        MaterialsProject ID number for the material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    k_dens_start : float, optional
        Starting kpoint density in 1/Angstrom. Defaults to 1.0
    k_dens_incr : float, optional
        Increment for the kpoint convergence. Can be set quite small since
        there is a check in place to see if a new mesh is actually constructed
        for each density. Defaults to 0.1.
    n_converge : int, optional
        Number of calculations that have to be inside the convergence
        threshold for convergence to be reached. Defaults to 3.
    high_level_db : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    """

    _fw_name = 'Start Encut Convergence'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file', 'k_dens_start', 'k_dens_incr', 'n_converge',
                       'high_level_db']

    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import converge_kpoints_swf
        mp_id = self.get('mp_id')
        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        k_dens_start = self.get('k_dens_start', 1.0)
        k_dens_incr = self.get('k_dens_incr', 0.1)
        n_converge = self.get('n_converge', 3)
        
        nav_structure = StructureNavigator(
            db_file=db_file, 
            high_level=hl_db)
        data = nav_structure.get_bulk_from_db(
            mp_id=mp_id, 
            functional=functional)
        
        stop_convergence = data.get('k_dens_info')
        
        if not stop_convergence:
            structure = Structure.from_dict(data.get('structure_equiVol'))
            comp_params = data.get('comp_parameters', {})
            SWF = converge_kpoints_swf(structure=structure,
                                       flag=mp_id,
                                       comp_parameters=comp_params,
                                       functional=functional,
                                       spec=fw_spec, 
                                       k_dens_start=k_dens_start,
                                       k_dens_incr=k_dens_incr, 
                                       n_converge=n_converge,
                                       db_file=db_file)

            return FWAction(detours=SWF, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)
        
@explicit_serialize
class FT_UpdateELists(FiretaskBase):
    """Fetch information about the last vasp calc and update the lists.
    
    Used with FT_KpointsConvo to converge k-mesh density using total energy
    calculations with vasp. This Firetasks reads the total energy from the
    StaticFW (atomate) with the task_label 'calc_name'. The shared date entry
    for the convergence is then identified via a tag and updated with the new
    total energy.
    
    Parameters
    ----------
    tag : str
        String from a uuid4 to identify the shared data entry in the database.
    calc_label : str
        A label for the specific calculation derived from the tag.
    db_file : str, optional
        Full path to the db.json file detailing access to the database.
        Defaults to '>>db_file<<' to use with env_chk.
    
    """
    
    _fw_name = 'Update Bulk Modulus Lists'
    required_params = ['tag', 'calc_label']
    optional_params = ['db_file']

    def run_task(self, fw_spec):
        
        tag = self.get('tag')
        calc_label = self.get('calc_label')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        nav = Navigator(db_file=db_file)
        # Get energy from last vasp run
        vasp_calc = nav.find_data(
            collection='tasks', 
            fltr={'task_label': calc_label})
        energy = vasp_calc['output']['energy']
        
        # Update data array in the database
        nav.update_data(
            collection='kpoints_data_sharing',
            fltr={'tag': tag},
            new_values={'$push': {'E_list': energy}})


@explicit_serialize
class FT_KpointsConvo(FiretaskBase):
    """Converge the kpoint density by monitoring the total energy.
    
    Runs a series of StaticFWs from atomate with increasing k-mesh density
    and monitors the resulting total energies until convergence is observed.
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-point grids.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
    tag : str
        String from a uuid4 to identify the shared data entry in the database
        and group the deformations calculations together.
    k_dens_start : int, optional
        Starting k-mesh density value for the first run. Defaults to 500.
    k_dens_incr : int, optional
        Increment for the k-mesh density during the convergence.
        Defaults to 50.
    n_converge : int, optional
        Number of calculations that have to show the same energy as the last
        one as to signify convergence, Defaults to 3.
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
        
    Returns
    -------
    FWActions that produce detour subworkflows until convergence is reached.
    """
    
    _fw_name = 'Kmesh density convergence'
    required_params = ['structure', 'comp_params', 'tag', 'flag',
                       'functional']
    optional_params = ['db_file', 'k_dens_start', 'k_dens_incr', 'n_converge',
                       'file_output', 'output_dir', 'remote_copy', 'server', 
                        'user', 'port', 'high_level_db']

    def run_task(self, fw_spec):
        
        n_converge = self.get('n_converge', 3)
        k_dens_start = self.get('k_dens_start', 1.0)
        k_dens_incr = self.get('k_dens_incr', 0.1)
        file_output = self.get('file_output', False)
        output_dir = self.get('output_dir', None)
        remote_copy = self.get('remote_copy', False)
        server = self.get('server', None)
        user = self.get('user', None)
        port = self.get('port', None)

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)
            
        struct = self['structure']
        comp_params = self['comp_params']
        tag = self['tag']
        flag = self['flag']
        functional = self['functional']
        
        # Get the relative energy tolerance (eV/atom) and the absolute tol.
        E_tolerance = comp_params.get('energy_tolerence', 0.001)
        E_tol = E_tolerance*struct.num_sites
        
        # Get the data arrays from the database (returns None when not there)
        nav = Navigator(db_file=db_file)
        data = nav.find_data(
            collection='kpoints_data_sharing',
            fltr={'tag': tag}
        )

        if data:
            E_list = data.get('E_list')
            k_dens_list = data.get('k_dens_list')
        else:
            E_list = None
            k_dens_list = None
        
        if E_list is None:
            
            label = tag+' calc 0'
            comp_params['k_dens'] = k_dens_start
            vis = get_custom_vasp_static_settings(struct, comp_params,
                                                  'bulk_from_scratch')
            kpoints = vis.kpoints
            k_dens_list = [k_dens_start]
                
            RunVASP_FW = StaticFW(structure=struct, vasp_input_set=vis,
                                  name=label)
            
            ParseAndUpdate_FW = Firework(
                [FT_UpdateELists(
                    tag=tag,
                    calc_label=label,
                    db_file=db_file),
                 FT_KpointsConvo(
                    structure=struct,
                    comp_params=comp_params,
                    tag=tag,
                    flag=flag,
                    functional=functional,
                    db_file=db_file,
                    k_dens_incr=k_dens_incr,
                    k_dens_start=k_dens_start,
                    n_converge=n_converge,
                    file_output=file_output,
                    output_dir=output_dir,
                    remote_copy=remote_copy,
                    server=server,
                    user=user,
                    port=port)],
                name='Update Energy Lists and Loop')
            
            K_convo_WF = Workflow([RunVASP_FW, ParseAndUpdate_FW],
                                  {RunVASP_FW: [ParseAndUpdate_FW]},
                                  name='Kpoint Convergence Loop')
            
            K_convo_WF = add_modify_incar(K_convo_WF)
            
            # Set up the entry for the data arrays in the database
            formula=struct.composition.reduced_formula
            set_data = {'tag': tag,
                        'chem_formula': formula,
                        'created_on': str(datetime.now()),
                        'k_dens_list': k_dens_list,
                        'last_mesh': kpoints.kpts,
                        'E_list': [],
                        'k-meshes': [kpoints.as_dict()]}
            nav.insert_data(
                collection='kpoints_data_sharing',
                data=set_data
            )
            
            return FWAction(detours=K_convo_WF)

        else:
            if is_list_converged(E_list, E_tol, n_converge):
                final_k_dens = k_dens_list[-n_converge]
                final_E = E_list[-n_converge]
                print('')
                print(' Convergence reached for total energy per atom.')
                print(' Final k_dens = {}; Final energy = {} eV;'
                      .format(final_k_dens, final_E))
                print('')
                print('')
        
                out_dict = {'k_dens_info': 
                                {'final_k_dens': final_k_dens,
                                 'final_energy': final_E,
                                 'energy_list': E_list,
                                 'k_dens_list': k_dens_list,
                                 'Energy_tol_abs': E_tol,
                                 'Energy_tol_rel': E_tolerance},
                            'total_energy@equiVol': final_E,
                            'comp_parameters.k_dens': final_k_dens}

                nav_high = Navigator(db_file=db_file, high_level=hl_db)
                nav_high.update_data(
                    collection=functional+'.bulk_data',
                    fltr={'mpid': flag},
                    new_values={'$set': out_dict},
                    upsert=True)

                nav.update_data(
                    collection='kpoints_data_sharing',
                    fltr={'tag': tag},
                    new_values={'$set': {'final_k_dens': final_k_dens,
                                         'final_energy': final_E,
                                         'energy_list': E_list,
                                         'k_dens_list': k_dens_list,
                                         'Energy_tol_abs': E_tol,
                                         'Energy_tol_rel': E_tolerance}})

                # Handle file output:
                if file_output:                 
                    write_FT = FileWriteTask(
                        files_to_write=[{'filename': flag+'_kpts_out.txt',
                                         'contents': pformat(out_dict)}])
                    copy_FT = copy_output_files(
                        file_list=[flag+'_kpts_out.txt'],
                        output_dir=output_dir,
                        remote_copy=remote_copy,
                        server=server,
                        user=server,
                        port=port)
                    FW = Firework([write_FT, copy_FT],
                                  name='Copy KpointsConvo SWF results')
                    WF = Workflow.from_Firework(
                        FW, 
                        name='Copy KpointsConve SWF results')

                    return FWAction(update_spec=fw_spec, detours=WF)
                else:  
                    return FWAction(update_spec=fw_spec)
            
            else:
                calc_nr = len(E_list)
                label = tag + ' calc '+str(calc_nr)
                
                last_mesh = data['last_mesh']
                k_dens = k_dens_list[-1] + k_dens_incr
                Mesh = MeshFromDensity(structure=struct,
                                       target_density=k_dens,
                                       compare_density=k_dens_list[-1])
                
                while Mesh.are_meshes_the_same():
                    k_dens = k_dens + k_dens_incr
                    Mesh = MeshFromDensity(structure=struct,
                                       target_density=k_dens,
                                       compare_density=k_dens_list[-1])
                comp_params['k_dens'] = k_dens
                vis = get_custom_vasp_static_settings(struct, comp_params,
                                                      'bulk_from_scratch')
                kpoints = vis.kpoints

                RunVASP_FW = StaticFW(structure=struct, vasp_input_set=vis,
                                      name=label)
            
                ParseAndUpdate_FW = Firework(
                    [FT_UpdateELists(
                         tag=tag,
                         calc_label=label,
                         db_file=db_file),
                     FT_KpointsConvo(
                        structure=struct,
                        comp_params=comp_params,
                        tag=tag,
                        flag=flag,
                        functional=functional,
                        db_file=db_file,
                        k_dens_incr=k_dens_incr,
                        k_dens_start=k_dens_start,
                        n_converge=n_converge,
                        file_output=file_output,
                        output_dir=output_dir,
                        remote_copy=remote_copy,
                        server=server,
                        user=user,
                        port=port)],
                    name='Update Energy Lists and Loop')

                K_convo_WF = Workflow([RunVASP_FW, ParseAndUpdate_FW],
                                      {RunVASP_FW: [ParseAndUpdate_FW]},
                                      name='Kpoint Convergence Loop')
                
                K_convo_WF = add_modify_incar(K_convo_WF)
                
                # Update Database entry for Encut list
                nav.update_data(
                    collection='kpoints_data_sharing',
                    fltr={'tag': tag},
                    new_values={'$push': {'k_dens_list': k_dens,
                                          'k-meshes': kpoints.as_dict()},
                                '$set': {'last_mesh': last_mesh}})

                return FWAction(detours=K_convo_WF)