"""Firetasks for converging the energy cutof in the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""
from datetime import datetime
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPStaticSet
from fireworks import FWAction, FiretaskBase, Firework, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.powerups import add_modify_incar
from triboflow.helper_functions import  GetCustomVaspStaticSettings, GetDB, \
    IsListConverged, GetBulkFromDB, GetHighLevelDB, GetGeneralizedKmesh


@explicit_serialize
class FT_StartKPointConvo(FiretaskBase):
    _fw_name = 'Start Encut Convergence'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file', 'k_dens_start', 'k_dens_incr', 'n_converge']
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import ConvergeKpoints_SWF
        mp_id = self.get('mp_id')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        k_dens_start = self.get('k_dens_start', 35)
        k_dens_incr = self.get('k_dens_incr', 0.5)
        n_converge = self.get('n_converge', 3)
        
        data = GetBulkFromDB(mp_id, db_file, functional)
        
        stop_convergence = data.get('kpoint_info')
        
        if not stop_convergence:
            structure = Structure.from_dict(data.get('structure_equiVol'))
            comp_params = data.get('comp_parameters', {})
            SWF = ConvergeKpoints_SWF(structure, comp_params,
                                    mp_id = mp_id, functional = functional,
                                    spec=fw_spec, 
                                    k_dens_start=k_dens_start,
                                    k_dens_incr=k_dens_incr, 
                                    n_converge=n_converge, db_file=db_file)
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
        DB = GetDB(db_file)
        
        #Get energy from last vasp run
        vasp_calc = DB.tasks.find_one({'task_label': calc_label})
        energy = vasp_calc['output']['energy']
        
        #update data array in the database
        DB.coll = DB['kpoints_data_sharing']
        DB.coll.update_one({'tag': tag},
                           {'$push': {'E_list': energy}})

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
    required_params = ['structure', 'comp_params', 'tag', 'mp_id',
                       'functional']
    optional_params = ['db_file', 'k_dens_start', 'k_dens_incr', 'n_converge']
    def run_task(self, fw_spec):
        
        n_converge = self.get('n_converge', 3)
        k_dens_start = self.get('k_dens_start', 35)
        k_dens_incr = self.get('k_dens_incr', 0.5)
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
            
        struct = self['structure']
        comp_params = self['comp_params']
        tag = self['tag']
        
        #get the relative energy tolerance (eV/atom) and the absolute tol.
        E_tolerance = comp_params.get('energy_tolerence', 0.001)
        E_tol = E_tolerance*struct.num_sites
        
        #get the data arrays from the database (returns None when not there)
        DB = GetDB(db_file)
        DB.coll = DB['kpoints_data_sharing']
        data = DB.coll.find_one({'tag': tag})
        if data:
            E_list = data.get('E_list')
            k_dens_list = data.get('k_dens_list')
        else:
            E_list = None
            k_dens_list = None
        
        if E_list is None:
            
            label = tag+' calc 0'
            vis, uis, vdw = GetCustomVaspStaticSettings(struct, comp_params,
                                                        'bulk_from_scratch')
            kpoints = Kpoints.automatic_gamma_density(struct, k_dens_start)
            #kpoints = GetGeneralizedKmesh(struct, k_dens_start)
            k_dens_list = [k_dens_start]
            vis = MPStaticSet(struct, user_incar_settings=uis,
                  user_kpoints_settings=kpoints)
                
            RunVASP_FW = StaticFW(structure=struct, vasp_input_set=vis,
                                  name=label)
            
            ParseAndUpdate_FW = Firework([FT_UpdateELists(tag = tag,
                                                          calc_label = label,
                                                          db_file = db_file),
                                          FT_KpointsConvo(structure = struct,
                                            comp_params = comp_params,
                                            tag = tag,
                                            mp_id = self['mp_id'],
                                            functional = self['functional'],
                                            db_file = db_file,
                                            k_dens_incr = k_dens_incr,
                                            k_dens_start = k_dens_start,
                                            n_converge = n_converge)],
                                         name='Update Energy Lists and Loop')
            
            K_convo_WF = Workflow([RunVASP_FW, ParseAndUpdate_FW],
                                  {RunVASP_FW: [ParseAndUpdate_FW]},
                                  name = 'Kpoint Convergence Loop')
            #Have to use this as a workaround as pymatgen forces ISMEAR=0
            #on a custom generated Kgrid, even if tetrahedra are provided.
            #K_convo_WF = add_modify_incar(K_convo_WF, {'incar_update': uis})
            
            #set up the entry for the data arrays in the database
            formula=struct.composition.reduced_formula
            set_data = {'tag': tag,
                        'chem_formula': formula,
                        'created_on': str(datetime.now()),
                        'k_dens_list': k_dens_list,
                        'last_mesh': kpoints.kpts,
                        'E_list': [],
                        'k-meshes': [kpoints.as_dict()]}
            DB.coll.insert_one(set_data)
            
            return FWAction(detours = K_convo_WF)
        
        else:
            if IsListConverged(E_list, E_tol, n_converge):
                final_k_dens = k_dens_list[-n_converge]
                final_E = E_list[-n_converge]
                print('')
                print(' Convergence reached for total energy per atom.')
                print(' Final k_dens = {}; Final energy = {}eV;'
                      .format(final_k_dens, final_E))
                print('')
                print('')
                mp_id = self.get('mp_id')
                functional = self.get('functional')
        
                tribo_db = GetHighLevelDB(db_file)
        
                hl_coll = tribo_db[functional+'.bulk_data']
                hl_coll.update_one({'mpid': mp_id},
                                   {'$set': {'k_dens_info': 
                                              {'final_k_dens': final_k_dens,
                                               'final_energy': final_E,
                                               'energy_list': E_list,
                                               'k_dens_list': k_dens_list,
                                               'Energy_tol_abs': E_tol,
                                               'Energy_tol_rel': E_tolerance},
                                    'total_energy@equiVol': final_E,
                                    'comp_parameters.k_dens': final_k_dens}})
                
                DB.coll.update_one({'tag': tag},
                               {'$set': {'final_k_dens': final_k_dens,
                                         'final_energy': final_E,
                                         'energy_list': E_list,
                                         'k_dens_list': k_dens_list,
                                         'Energy_tol_abs': E_tol,
                                         'Energy_tol_rel': E_tolerance}})
                
                return FWAction(update_spec = fw_spec)
            
            else:
                calc_nr = len(E_list)
                label = tag+' calc '+str(calc_nr)
                vis, uis, vdw = GetCustomVaspStaticSettings(struct, comp_params,
                                                        'bulk_from_scratch')
                last_mesh = data['last_mesh']
                k_dens = k_dens_list[-1]+k_dens_incr
                kpoints = Kpoints.automatic_gamma_density(struct, k_dens)
                #kpoints = GetGeneralizedKmesh(struct, k_dens)
                while kpoints.kpts == last_mesh:
                    k_dens = k_dens + k_dens_incr
                    kpoints = Kpoints.automatic_gamma_density(struct, k_dens)
                    #kpoints = GetGeneralizedKmesh(struct, k_dens)
                vis = MPStaticSet(struct, user_incar_settings=uis,
                                  user_kpoints_settings=kpoints)
                
                print(comp_params)
                print(vis.incar)

                RunVASP_FW = StaticFW(structure=struct, vasp_input_set=vis,
                                      name=label)
            
                ParseAndUpdate_FW = Firework([FT_UpdateELists(tag = tag,
                                                          calc_label = label,
                                                          db_file = db_file),
                                          FT_KpointsConvo(structure = struct,
                                            comp_params = comp_params,
                                            tag = tag,
                                            mp_id = self['mp_id'],
                                            functional = self['functional'],
                                            db_file = db_file,
                                            k_dens_incr = k_dens_incr,
                                            k_dens_start = k_dens_start,
                                            n_converge = n_converge)],
                                         name='Update Energy Lists and Loop')
            
                K_convo_WF = Workflow([RunVASP_FW, ParseAndUpdate_FW],
                                      {RunVASP_FW: [ParseAndUpdate_FW]},
                                      name = 'Kpoint Convergence Loop')
                #Have to use this as a workaround as pymatgen forces ISMEAR=0
                #on a custom generated Kgrid, even if tetrahedra are provided.
                K_convo_WF = add_modify_incar(K_convo_WF, {'incar_update': uis})
                
                #Update Database entry for Encut list
                DB.coll.update_one({'tag': tag},
                               {'$push': {'k_dens_list': k_dens,
                                          'k-meshes': kpoints.as_dict()},
                                '$set': {'last_mesh': kpoints.kpts}})
                return FWAction(detours=K_convo_WF)

