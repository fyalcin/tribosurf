"""Firetasks for converging the energy cutof in the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from datetime import datetime
from pprint import pprint, pformat

import numpy as np
import pymongo
from pymatgen.core.structure import Structure
from fireworks import FWAction, FiretaskBase, Firework, Workflow, FileWriteTask
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.workflows.base.bulk_modulus import get_wf_bulk_modulus

from triboflow.utils.vasp_tools import get_custom_vasp_static_settings, get_emin
from triboflow.utils.check_convergence import is_list_converged
from triboflow.utils.file_manipulation import copy_output_files
from triboflow.utils.database import Navigator



@explicit_serialize
class FT_StartEncutConvo(FiretaskBase):
    _fw_name = 'Start Encut Convergence'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file']
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import ConvergeEncut_SWF
        mp_id = self.get('mp_id')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        nav = Navigator(db_file)
        data = nav.find_data(functional+'.bulk_data', {'mpid': mp_id})
        
        stop_convergence = data.get('encut_info')
        
        if not stop_convergence:
            structure = Structure.from_dict(data.get('structure_fromMP'))
            comp_params = data.get('comp_parameters', {})
            SWF = ConvergeEncut_SWF(structure = structure, flag = mp_id,
                                    comp_parameters = comp_params,
                                    functional = functional,
                                    print_help = False)
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
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        nav = Navigator(db_file)
        results = nav.find_many_data(nav.db.eos, {'formula_pretty': formula})

        # Get the first element of the ordered results
        ordered_results = results.sort('created_at', pymongo.DESCENDING)[0]

        BM = ordered_results['bulk_modulus']
        V0 = ordered_results['results']['v0']
        
        # Update data arrays in the database
        nav.update_data('BM_data_sharing', {'tag': tag}, 
                        {'$push': {'BM_list': BM, 'V0_list': V0}})


@explicit_serialize
class FT_EnergyCutoffConvo(FiretaskBase):
    """Converge the encut for a material via fits to an EOS, raising encut.
    
    Uses the get_bulk_modulus workflow of atomate to fit Birch-Murnaghen EOS
    for increasing values of the d = mmdb.collection.find_one({"task_label": {"$regex": "{} bulk_modulus*".format(tag)}})energy cutoff. Once bulk modulus and
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
        Starting encut value for the first run. Defaults to the largest EMIN
        in the POTCAR.
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
    required_params = ['structure', 'comp_params', 'tag', 'flag',
                       'functional']
    optional_params = ['deformations', 'n_converge', 'encut_start',
                       'encut_incr', 'db_file', 'file_output', 'output_dir',
                       'remote_copy', 'server', 'user', 'port']
    def run_task(self, fw_spec):
        deforms = []
        for i in np.arange(0.9, 1.1, 0.05):
            dm=np.eye(3)*i
            deforms.append(dm)  
        n_converge = self.get('n_converge', 3)
        encut_start = self.get('encut_start', None)
        encut_incr = self.get('encut_incr', 25)
        deformations = self.get('deformations')
        file_output = self.get('file_output', False)
        output_dir = self.get('output_dir', None)
        remote_copy = self.get('remote_copy', False)
        server = self.get('server', None)
        user = self.get('user', None)
        port = self.get('port', None)
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        if not deformations:
            deformations = deforms
            
        struct = self['structure']
        comp_params = self['comp_params']
        tag = self['tag']
        
        V0_tolerance = comp_params.get('volume_tolerence', 0.001)
        BM_tolerance = comp_params.get('BM_tolerence', 0.01)
        #uks = {'reciprocal_density': 1000}
        
        
        # Get the data arrays from the database (returns None when not there)
        nav = Navigator(db_file)
        data = nav.find_data('BM_data_sharing', {'tag': tag})
        
        if data:
            BM_list = data.get('BM_list')
            V0_list = data.get('V0_list')
            Encut_list = data.get('Encut_list')
        else:
            BM_list = None
            V0_list = None
            Encut_list = None
        
        if BM_list is None:
            vis = get_custom_vasp_static_settings(struct, comp_params,
                                              'bulk_from_scratch')            
            if not encut_start:
                # Get the largest EMIN value of the potcar and round up to the
                # next whole 25.
                emin = get_emin(vis.potcar)
                encut_start = int(25 * np.ceil(emin/25))
            
            vis.user_incar_settings.update({'ENCUT': encut_start})
            Encut_list = [encut_start]

            BM_WF = get_wf_bulk_modulus(struct, deformations,
                                        vasp_input_set=vis,
                                        vasp_cmd=VASP_CMD, db_file=db_file,
                                        #user_kpoints_settings=uks,
                                        eos='birch_murnaghan', tag=tag)
            
            formula=struct.composition.reduced_formula
            UAL_FW = Firework([FT_UpdateBMLists(formula=formula, tag=tag),
                               FT_EnergyCutoffConvo(structure = struct,
                                         comp_params = comp_params,
                                         tag = tag,
                                         flag = self['flag'],
                                         functional = self['functional'],
                                         db_file = db_file,
                                         encut_incr = encut_incr,
                                         encut_start = encut_start,
                                         file_output = file_output,
                                         output_dir = output_dir,
                                         remote_copy = remote_copy,
                                         server = server,
                                         user = user,
                                         port = port)],
                              name='Update BM Lists and Loop')
            
            BM_WF.append_wf(Workflow.from_Firework(UAL_FW), BM_WF.leaf_fw_ids)
            # Use add_modify_incar powerup to add KPAR and NCORE settings
            # based on env_chk in my_fworker.yaml
            BM_WF = add_modify_incar(BM_WF)
            
            # Set up the entry for the data arrays in the database
            nav.insert_data('BM_data_sharing', 
                            {'tag': tag,
                             'chem_formula': formula,
                             'created_on': str(datetime.now()),
                             'Encut_list': Encut_list,
                             'BM_list': [],
                             'V0_list': []})
            
            return FWAction(detours=BM_WF)
        
        else:
            BM_tol = BM_list[-1]*BM_tolerance
            V0_tol = V0_list[-1]*V0_tolerance
            if (is_list_converged(BM_list, BM_tol, n_converge)
            and is_list_converged(V0_list, V0_tol, n_converge)):
                final_encut = Encut_list[-n_converge]
                final_BM = BM_list[-n_converge]
                final_V0 = V0_list[-n_converge]
                flag = self.get('flag')
                functional = self.get('functional')
        
                scaled_structure = struct.copy()
                scaled_structure.scale_lattice(final_V0)
                struct_dict = scaled_structure.as_dict()
                
                print('')
                print(' Convergence reached for BM and cell volume.')
                print(' Final encut = {} eV; Final BM = {} GPa; Final Volume = {} Angstrom³'
                      .format(final_encut, final_BM, final_V0))
                print('')
                print(' The scaled output structure is:\n')
                pprint(struct_dict)
        
        
                output_dict = {'encut_info': 
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
                               'equilibrium_volume': final_V0,
                               'bulk_moduls': final_BM,
                               'comp_parameters.encut': final_encut,
                               'structure_equiVol': struct_dict}
                
                nav_high = Navigator(db_file, high_level='triboflow')
                nav_high.update_data(
                    functional+'.bulk_data',
                    {'mpid': flag},
                    {'$set': output_dict},
                    upsert=True)

                nav.update_data(
                    'BM_data_sharing',
                    {'mpid': tag},
                    {'$set': {'final_encut': final_encut,
                              'final_BM': final_BM,
                              'final_volume': final_V0,
                              'BM_tol_abs': BM_tol,
                              'BM_tol_rel': BM_tolerance,
                              'V0_tol_abs': V0_tol,
                              'V0_tol_rel': V0_tolerance}})

                # handle file output:
                if file_output:                 
                    write_FT = FileWriteTask(
                        files_to_write=[{'filename': flag+'_output_dict.txt',
                                         'contents': pformat(output_dict)}])

                    copy_FT = copy_output_files(
                        file_list = [flag+'_output_dict.txt'],
                        output_dir = output_dir,
                        remote_copy = remote_copy,
                        server = server,
                        user = server,
                        port = port)

                    FW = Firework(
                        [write_FT, copy_FT],
                        name = 'Copy Encut SWF results')

                    WF = Workflow.from_Firework(
                        FW,
                        name = 'Copy Encut SWF results')

                    return FWAction(update_spec = fw_spec, detours = WF)
                else:  
                    return FWAction(update_spec = fw_spec)
            
            vis = get_custom_vasp_static_settings(struct, comp_params,
                                              'bulk_from_scratch')
            encut = Encut_list[-1]+encut_incr
            vis.user_incar_settings.update({'ENCUT': encut})
            
            BM_WF = get_wf_bulk_modulus(struct, deformations,
                                        vasp_input_set=vis,
                                        vasp_cmd=VASP_CMD, db_file=DB_FILE,
                                        #user_kpoints_settings=uks,
                                        eos='birch_murnaghan', tag=tag)
            
            formula=struct.composition.reduced_formula
            UAL_FW = Firework(
                [FT_UpdateBMLists(formula=formula, tag=tag),
                FT_EnergyCutoffConvo(structure=struct,
                                     comp_params=comp_params,
                                     tag=tag,
                                     flag=self['flag'],
                                     functional=self['functional'],
                                     db_file=db_file,
                                     encut_incr=encut_incr,
                                     encut_start=encut_start,
                                     file_output=file_output,
                                     output_dir=output_dir,
                                     remote_copy=remote_copy,
                                     server=server,
                                     user=user,
                                     port=port)],
                name='Update BM Lists and Loop')
            
            BM_WF.append_wf(Workflow.from_Firework(UAL_FW), BM_WF.leaf_fw_ids)
            # Use add_modify_incar powerup to add KPAR and NCORE settings
            # based on env_chk in my_fworker.yaml
            BM_WF = add_modify_incar(BM_WF)
            
            # Update Database entry for Encut list
            nav.update_data(
                'BM_data_sharing',
                {'tag': tag}, 
                {'$push': {'Encut_list': encut}})

            return FWAction(detours=BM_WF)

