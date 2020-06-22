"""Firetasks for converging the energy cutof in the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""
import numpy as np
from datetime import datetime
from pymatgen.core.structure import Structure
from fireworks import FWAction, FiretaskBase, Firework, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.workflows.base.bulk_modulus import get_wf_bulk_modulus
from triboflow.helper_functions import GetLastBMDatafromDB, \
    GetCustomVaspStaticSettings, GetDB, IsListConverged, GetBulkFromDB, \
    GetHighLevelDB


@explicit_serialize
class FT_StartEncutConvo(FiretaskBase):
    _fw_name = 'Start Encut Convergence'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file', 'deformations', 'encut_start', 'encut_incr',
                       'n_converge']
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import ConvergeEncut_SWF
        mp_id = self.get('mp_id')
        functional = self.get('functional')
        db_file = self.get('db_file', env_chk('>>db_file<<', fw_spec))
        deformations = self.get('deformations')
        encut_start = self.get('encut_start', 200)
        encut_incr = self.get('encut_incr', 25)
        n_converge = self.get('n_converge', 3)
        
        data = GetBulkFromDB(mp_id, db_file, functional)
        
        stop_convergence = data.get('encut_info')
        
        if not stop_convergence:
            structure = Structure.from_dict(data.get('structure_fromMP'))
            comp_params = data.get('comp_parameters', {})
            SWF = ConvergeEncut_SWF(structure, comp_params,
                                    mp_id = mp_id, functional = functional,
                                    spec=fw_spec, 
                                    deformations=deformations,
                                    encut_start=encut_start,
                                    encut_incr=encut_incr, 
                                    n_converge=n_converge, db_file=db_file)
            return FWAction(detours=SWF, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)
        

@explicit_serialize
class FT_UpdateBMLists(FiretaskBase):
    """Fetch information about the EOS fit from the DB and update the lists.
    
    Used with FT_EnergyCutoffConvo to converge energy cutoffs using the
    get_wf_bulk_modulus workflow of atomate. This Firetasks reads the
    neccessary information from the eos collection of the database. Since no
    usefull tag is placed, the identification of the correct entry is done
    by the chemical formula and the timestamp. The shared date entry for the
    convergence is then identified via a tag and updated with the new
    equilibrium volume and the bulk modulus.
    
    Parameters
    ----------
    formula : str
        Chemical formula on the material to be matched with the database.
    tag : str
        String from a uuid4 to identify the shared data entry in the database.
    db_file : str, optional
        Full path to the db.json file detailing access to the database.
        Defaults to '>>db_file<<' to use with env_chk.
    
    """
    
    _fw_name = 'Update Bulk Modulus Lists'
    required_params = ['formula', 'tag']
    optional_params = ['db_file']
    def run_task(self, fw_spec):
        formula = self.get('formula')
        tag = self.get('tag')
        db_file = self.get('db_file', env_chk('>>db_file<<', fw_spec))
        
        results = GetLastBMDatafromDB(formula, db_file)
        
        BM = results['bulk_modulus']
        V0 = results['results']['v0']
        
        #update data arrays in the database
        DB = GetDB(db_file)
        DB.coll = DB['BM_data_sharing']
        DB.coll.update_one({'tag': tag},
                           {'$push': {'BM_list': BM,
                                      'V0_list': V0}})


@explicit_serialize
class FT_EnergyCutoffConvo(FiretaskBase):
    """Converge the encut for a material via fits to an EOS, raising encut.
    
    Uses the get_bulk_modulus workflow of atomate to fit Birch-Murnaghen EOS
    for increasing values of the energy cutoff. Once bulk modulus and
    equilibrium volume are converged, the subsequent detours are stopped and
    the convergence data is passed on.
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-pint grids.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
    tag : str
        String from a uuid4 to identify the shared data entry in the database
        and group the deformations calculations together.
    deformations: list of lists, optional
        List of deformation matrices for the fit to the EOS. Defaults to None,
        which results in 5 volumes from 90% to 110% of the initial volume.
    encut_start : float, optional
        Starting encut value for the first run. Defaults to 200.
    encut_incr : float, optional
        Increment for the encut during the convergence. Defaults to 25.
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
    
    _fw_name = 'Energy Cutoff Convergence'
    required_params = ['structure', 'comp_params', 'tag', 'mp_id',
                       'functional']
    optional_params = ['deformations', 'n_converge', 'encut_start',
                       'encut_incr', 'db_file']
    def run_task(self, fw_spec):
        deforms = []
        for i in np.arange(0.9, 1.1, 0.05):
            dm=np.eye(3)*i
            deforms.append(dm)  
        n_converge = self.get('n_converge', 3)
        encut_start = self.get('encut_start', 200)
        encut_incr = self.get('encut_incr', 25)
        deformations = self.get('deformations', deforms)
        db_file = self.get('db_file', env_chk('>>db_file<<', fw_spec))
            
        struct = self['structure']
        comp_params = self['comp_params']
        tag = self['tag']
        
        V0_tolerance = comp_params.get('volume_tolerence', 0.001)
        BM_tolerance = comp_params.get('BM_tolerence', 0.01)
        
        
        #get the data arrays from the database (returns None when not there)
        DB = GetDB(db_file)
        DB.coll = DB['BM_data_sharing']
        data = DB.coll.find_one({'tag': tag})
        if data:
            BM_list = data.get('BM_list')
            V0_list = data.get('V0_list')
            Encut_list = data.get('Encut_list')
        else:
            BM_list = None
            V0_list = None
            Encut_list = None
        
        if BM_list is None:
            vis, uis, vdw = GetCustomVaspStaticSettings(struct, comp_params,
                                                        'bulk_from_scratch')
            uis['ENCUT'] = encut_start
            Encut_list = [encut_start]
            BM_WF = get_wf_bulk_modulus(struct, deformations,
                                        vasp_input_set=None,
                                        vasp_cmd=VASP_CMD, db_file=db_file,
                                        user_kpoints_settings=None,
                                        eos='birch_murnaghan', tag=tag,
                                        user_incar_settings=uis)
            
            formula=struct.composition.reduced_formula
            UAL_FW = Firework([FT_UpdateBMLists(formula=formula, tag=tag),
                               FT_EnergyCutoffConvo(structure = struct,
                                         comp_params = comp_params,
                                         tag = tag,
                                         mp_id = self['mp_id'],
                                         functional = self['functional'],
                                         db_file = db_file,
                                         encut_incr = encut_incr,
                                         encut_start = encut_start)],
                              name='Update BM Lists and Loop')
            
            BM_WF.append_wf(Workflow.from_Firework(UAL_FW), BM_WF.leaf_fw_ids)
            
            #set up the entry for the data arrays in the database
            set_data = {'tag': tag,
                        'chem_formula': formula,
                        'created_on': str(datetime.now()),
                        'Encut_list': Encut_list,
                        'BM_list': [],
                        'V0_list': []}
            DB.coll.insert_one(set_data)
            
            return FWAction(detours=BM_WF)
        
        else:
            BM_tol = BM_list[-1]*BM_tolerance
            V0_tol = V0_list[-1]*V0_tolerance
            if (IsListConverged(BM_list, BM_tol, n_converge)
            and IsListConverged(V0_list, V0_tol, n_converge)):
                final_encut = Encut_list[-n_converge]
                final_BM = BM_list[-n_converge]
                final_V0 = V0_list[-n_converge]
                print('')
                print(' Convergence reached for BM and cell volume.')
                print(' Final encut = {}eV; Final BM = {}GPa; Final Volume = {}'
                      .format(final_encut, final_BM, final_V0))
                print('')
                print('')
                mp_id = self.get('mp_id')
                functional = self.get('functional')
        
                scaled_structure = struct.copy()
                scaled_structure.scale_lattice(final_V0)
                struct_dict = scaled_structure.as_dict()
                struct_name = 'structure_equiVol'
        
                tribo_db = GetHighLevelDB(db_file)
        
                coll = tribo_db[functional+'.bulk_data']
                coll.update_one({'mpid': mp_id},
                                {'$set': {'encut_info': 
                                              {'final_encut': final_encut,
                                               'final_BM': final_BM,
                                               'final_volume': final_V0,
                                               'BM_list': BM_list,
                                               'V0_list': V0_list,
                                               'Encut_list': Encut_list,
                                               'BM_tol_abs': BM_tol,
                                               'BM_tol_rel': BM_tolerance,
                                               'V0_tol_abs': V0_tol,
                                               'V0_tol_rel': V0_tolerance},
                                        'encut': final_encut,
                                        'equilibrium_volume': final_V0,
                                        'bulk_moduls': final_BM,
                                        struct_name: struct_dict}})
                
                DB.coll.update_one({'tag': tag},
                               {'$set': {'final_encut': final_encut,
                                         'final_BM': final_BM,
                                         'final_volume': final_V0,
                                         'BM_tol_abs': BM_tol,
                                         'BM_tol_rel': BM_tolerance,
                                         'V0_tol_abs': V0_tol,
                                         'V0_tol_rel': V0_tolerance}})
                
                return FWAction(update_spec = fw_spec)
            
            else:
                vis, uis, vdw = GetCustomVaspStaticSettings(struct, comp_params,
                                                        'bulk_from_scratch')
            encut = Encut_list[-1]+encut_incr
            uis['ENCUT'] = encut
            
            BM_WF = get_wf_bulk_modulus(struct, deformations,
                                        vasp_input_set=None,
                                        vasp_cmd=VASP_CMD, db_file=DB_FILE,
                                        user_kpoints_settings=None,
                                        eos='birch_murnaghan', tag=tag,
                                        user_incar_settings=uis)
            
            formula=struct.composition.reduced_formula
            UAL_FW = Firework([FT_UpdateBMLists(formula=formula, tag=tag),
                               FT_EnergyCutoffConvo(structure = struct,
                                         comp_params = comp_params,
                                         tag = tag,
                                         mp_id = self['mp_id'],
                                         functional = self['functional'],
                                         db_file = db_file,
                                         encut_incr = encut_incr,
                                         encut_start = encut_start)],
                              name='Update BM Lists and Loop')
            
            BM_WF.append_wf(Workflow.from_Firework(UAL_FW), BM_WF.leaf_fw_ids)
            
            #Update Database entry for Encut list
            DB.coll.update_one({'tag': tag},
                               {'$push': {'Encut_list': encut}})
            return FWAction(detours=BM_WF)

