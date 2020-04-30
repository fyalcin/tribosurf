#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python
# -*- coding: utf-8 -*-
"""A collection of Firetasks to be used for the FireFlow project.

Created on Mon Mar 23 14:33:47 2020
@author: mwo
"""
from fireworks.core.firework import FWAction, FiretaskBase, Firework
from fireworks.utilities.fw_utilities import explicit_serialize

@explicit_serialize
class FT_LoopKpoints(FiretaskBase):
    _fw_name = 'Converge Kpoints'
    required_params = ['structure_loc', 'kinfo_loc', 'out_loc', 'comp_params',
                       'k_dist_start', 'k_dist_incr']
    optional_params = ['n_converge']
    def run_task(self, fw_spec):
        from CommonWorkflows import GetEnergy_SWF
        from HelperFunctions import GetValueFromNestedDict
   
        if 'n_converge' in self:
            n_converge = self['n_converge']
        else:
            n_converge = 3
            
        final_k_dist = GetValueFromNestedDict(fw_spec,
                                        self['kinfo_loc']+['final_k_distance'])
        if final_k_dist:
            k_str = '->'.join(self['out_loc']+['final_k_distance'])
            return FWAction(mod_spec=[{'_set': {k_str: final_k_dist}}])
        else:
            FW = Firework([FT_ConvergeKpoints, FT_LoopKpoints])
            return FWAction(detours=Workflow(FW))
            
@explicit_serialize
class FT_ConvergeKpoints(FiretaskBase):
    _fw_name = 'Converge Kpoints'
    required_params = ['structure_loc', 'kinfo_loc', 'out_loc', 'comp_params',
                       'k_dist_start', 'k_dist_incr', 'n_converge']
    def run_task(self, fw_spec):
        from CommonWorkflows import GetEnergy_SWF
        from HelperFunctions import GetValueFromNestedDict, IsEnergyConverged
        
        energies_loc = self['kinfo_loc']+['energy_array']
        energies = GetValueFromNestedDict(fw_spec, energies_loc)
        n_converge = self['n_converge']
        struct = GetValueFromNestedDict(fw_spec, self ['structure'])
        comp_parameters = GetValueFromNestedDict(fw_spec, self['comp_params'])
        #TODO: Decide the static_type based on the structure (bulk or slab)!
        static_type = 'bulk_from_scratch'
        out_loc = self['kinfo_loc']+['last_calc']
        to_pass = [self['kinfo_loc'][0]]
        spec = fw_spec
        etol = comp_parameters['energy_tolerance']*struct.num_sites
        k_dist = GetValueFromNestedDict(fw_spec, 
                                        self['kinfo_loc']+['k_distences'])[-1]
        
        if energies is None:
            k_dist = self['k_dist_start']
            comp_parameters['k_distance'] = k_dist
            k_str = '->'.join(self['kinfo_loc']+['k_distences'])
            return FWAction(detours=GetEnergy_SWF(struct, comp_parameters,
                                                  static_type, out_loc,
                                                  to_pass, spec,
                                                  push_energy_loc=energies_loc
                                                  ),
                            mod_spec=[{'_push': {k_str: k_dist}}])
        elif len(energies) <= n_converge:
            k_dist = k_dist + self['k_dist_incr']
            comp_parameters['k_distance'] = k_dist
            k_str = '->'.join(self['kinfo_loc']+['k_distences'])
            return FWAction(detours=GetEnergy_SWF(struct, comp_parameters,
                                                  static_type, out_loc,
                                                  to_pass, spec,
                                                  push_energy_loc=energies_loc
                                                  ),
                            mod_spec=[{'_push': {k_str: k_dist}}])
        else:
            if IsEnergyConverged(energies, etol, n_converge):
                final_k_dist = GetValueFromNestedDict(fw_spec, 
                            self['kinfo_loc']+['k_distences'])[-n_converge]
                k_str = '->'.join(self['kinfo_loc']+['final_k_distance'])
                return FWAction(mod_spec=[{'_set': {k_str: final_k_dist}}])
            
            
        
        

@explicit_serialize
class FT_GetEnergyFromDB(FiretaskBase):
    """ Make a database query to get the final energy for a VASP calculation.
    
    This might not be unambigous if you have more than one calculation with
    the same formula and task_label. The energy will be the one of the first
    matching calculation found in the database.
    
    Parameters
    ----------
    task_label : str
        Specifies what type of vasp calculation you are querying. E.g. 'static'
        or 'nscf uniform'.
    formula : str
        Reduced chemical formula for the compound you are looking for.
    out_loc : list of str
        Location in the spec as specified by a list of keys into which the
        energy should be written.
        
    Returns
    -------
    The final energy of a VASP calculation at a specified location in the spec.
    """
    _fw_name = 'Get Energy'
    required_params = ['label', 'formula', 'out_loc']
    def run_task(self, fw_spec):
        from atomate.vasp.database import VaspCalcDb
        db_file = fw_spec['_fw_env']['db_file']
        atomate_db = VaspCalcDb.from_db_file(db_file=db_file)
        
        label = self['label']
        formula = self['formula']
        out_loc = self['out_loc']
        
        energy = atomate_db.collection.find_one({'task_label': label,
                                                 'formula_pretty': formula})[
                                                     'output']['energy']
        out_str = '->'.join(out_loc)
        return FWAction(mod_spec=[{'_set': {out_str: energy}}])
    
@explicit_serialize
class FT_StartRelaxSubWorkflow(FiretaskBase):
    """Start a RelaxSubWorkflow as a detour.
    
    This Firetask is intendet to be used within a StartDetourWF_FW Firework.
    A Relax_SWF is started as a detour using a FWAction. The current fw_spec
    is passed on to the subworkflow.
    
    Parameters
    ----------
    structure_loc : list of str
        List of keys that point to the structure to be relaxed in the fw_spec.
    comp_parameters_loc : list of str
        List of keys that point to the computational parameters to be used for
        the relaxation in the fw_spec.
    relax_type : str
        specifies which relaxation should be performed.
    out_loc : list of str
        List of keys that point to the relaxed structure in the fw_spec.
    to_pass : list of str
        List of keys that each represent a location in the first level of the
        fw_spec and signifies which of those are to be passed to the next
        Firework. E.g. if the fw_sec = {'a': a, 'b': {'b1'; b1, 'b2': b2},
        'c': c}} and to_pass = ['a', 'b'], the spec in the next FW will be:
        {'a': a, 'b': {'b1'; b1, 'b2': b2}}
    
    See Also
    --------
    StartDetourWF_FW in CommonFireworks and Relax_SWF in CommonWorkflows.
    
    Returns
    -------
    FWAction with a detour starting a Relax_SWF subworkflow.
    """
    
    _fw_name = 'Starting Relaxation Subworkflow'
    required_params = ['structure_loc', 'comp_parameters_loc', 'relax_type',
                       'out_loc', 'to_pass']
    def run_task(self, fw_spec):
        from CommonWorkflows import Relax_SWF
        from HelperFunctions import GetValueFromNestedDict
        
        relax_type = self['relax_type']
        structure = GetValueFromNestedDict(fw_spec, self['structure_loc'])
        params = GetValueFromNestedDict(fw_spec, self['comp_parameters_loc'])
        out_loc = self['out_loc']
        to_pass = self['to_pass']
        
        SWF = Relax_SWF(structure, params, relax_type, out_loc, to_pass,
                        fw_spec)
        return FWAction(detours=SWF)
        
@explicit_serialize
class FT_ParseVaspOutput(FiretaskBase):
    """Parse VASP OUTCAR output and save it in a given location in the fw_spec.
    
    The OUTCAR is parsed using pymatgen Outcar object. The location of the
    vasp calculation folder can be directly specified, but probably the best
    way is to first copy the output files using atomate PassCalcLocs features.
    FireWorks _pass_job_info functionality is also supported.
    
    Parameters
    ----------
    out_loc : list of str
        List of keys that point to the location in the fw_spec that should
        contain the results that are parsed.
    job_dir : str, optional
        If job_dir is passed it is used to specify the directory where the
        VASP job was run. 
        
    Returns
    -------
    An output directory with OUTCAR results is set in the fw_spec using a
    FWAction.
    """
    _fw_name = 'Parse Vasp Output'
    required_params = ['out_loc']
    optional_params = ['job_dir', 'push_energy_loc']
    def run_task(self, fw_spec):
        import os.path
        from pymatgen.io.vasp.outputs import Outcar
        from pymatgen.core import Structure
        from HelperFunctions import SetInSpecFWAction
        
        if 'job_dir' in self:
            job_dir = self['job_dir']
        elif '_job_info' in fw_spec:
            job_dir = fw_spec['_job_info'][-1]['launch_dir']
        else:
            job_dir = '.'
            
        if os.path.isfile(job_dir+'/OUTCAR.relax2.gz'):
            outcar = Outcar(job_dir+'/OUTCAR.relax2.gz')
        elif os.path.isfile(job_dir+'/OUTCAR.gz'):
            outcar = Outcar(job_dir+'/OUTCAR.gz')
        elif os.path.isfile(job_dir+'/OUTCAR.relax2'):
            outcar = Outcar(job_dir+'/OUTCAR.relax2')
        elif os.path.isfile(job_dir+'/OUTCAR'):
            outcar = Outcar(job_dir+'/OUTCAR')
        else:
            raise IOError('OUTCAR could not be located in folder:\n{}'.format(
                job_dir))
        
        output_dict = {}
        output_dict['energy'] = outcar.final_energy
        output_dict['fermi_energy'] = outcar.efermi
        output_dict['local_charges'] = outcar.charge
        output_dict['local_magnetization'] = outcar.magnetization
        
        if self.get('push_energy_loc'):
            out_str = '->'.join(self['out_loc'])
            push_str = '->'.join(self['push_energy_loc'])
            FWA = FWAction(mod_spec=[{'_set': {out_str: output_dict}},
                                     {'_push': {push_str: outcar.final_energy}}
                                    ])
        else:
            FWA = SetInSpecFWAction(self['out_loc'], output_dict)
        return FWA
        
@explicit_serialize
class FT_StructFromVaspOutput(FiretaskBase):
    """Make a pymatgen Structure object out of VASP output files.
    
    Parameters
    ----------
    out_struct_loc: list of str
        location of the structure in the spec given by a list of keys.
    job_dir: str, optional
        directory where the VASP files (CONTCAR.relax.2.gz and
        OUTCAR.relax.2.gz) are located. This assumes that a double-relaxation
        run has been done before with custodian. If not defined, the files
        will be loaded from the current directory. Since the files might be on
        a different machine, better use CopyFilesFromCalcLoc before.
    ----------
    
    Returns
    -------
    structure : pymatgen.core.structure.Structure
        The structure is added to the fw_spec.
    -------
    
    Notes
    -----
    TODO: Maybe this should be generalized to except specific file names.
    However, this have to be optional parameters and should default to the
    current values of CONTCAR.relax.2.gz and OUTCAR.relax.2.gz!
    
    """
    
    _fw_name = 'Make Structure from Vasp Output'
    required_params = ['out_struct_loc']
    optional_params = ['job_dir']
    def run_task(self, fw_spec):
        from pymatgen.io.vasp.outputs import Outcar
        from pymatgen.core import Structure
        from HelperFunctions import WriteNestedDictFromList, UpdateNestedDict
        
        if 'job_dir' in self:
            job_dir = self['job_dir']
        elif '_job_info' in fw_spec:
            job_dir = fw_spec['_job_info'][-1]['launch_dir']
        else:
            job_dir = '.'
        outcar = Outcar(job_dir+'/OUTCAR.relax2.gz')
        structure = Structure.from_file(job_dir+'/CONTCAR.relax2.gz')
        
        mag_tuple = outcar.magnetization
        
        if mag_tuple:
            for i, site in enumerate(structure.sites):
                structure.replace(i, structure.species[i],
                                  properties={'magmom': mag_tuple[i]['tot']})
        
        output = WriteNestedDictFromList(self['out_struct_loc'],
                                         structure, d={})

        spec = fw_spec
        updated_spec = UpdateNestedDict(spec, output)
        return FWAction(update_spec = updated_spec)
    
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
    
    _fw_name = 'Add selctive dynamics'
    required_params = ['structure']
    optional_params = ['selective_dynamics_array']

    def run_task(self, fw_spec):
        import numpy as np
        from pymatgen.io.vasp.inputs import Poscar
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
        

@explicit_serialize
class FT_SpawnOptimizeFW(FiretaskBase):
    """Spawns a new optimization workflow using a structure saved in the spec.
    
    The structure has to be saved in the 'structures' dictionary within the
    fw_spec and has to have the 'to_optimize' key.
    """
    
    _fw_name = 'Spawn Optimize_Workflow'

    def run_task(self, fw_spec):
        from atomate.vasp.fireworks.core import OptimizeFW
        from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet
        
        structure = fw_spec['structures']['to_optimize']
        
        if fw_spec['Input_test.in']['functional'] == 'SCAN':
            vis = MPScanRelaxSet(structure)
        else:
            vis = MPRelaxSet(structure)
        
        #SCAN inut currently does not work for OptimizeFW!
        #new_optimize_wf = OptimizeFW(structure, vasp_input_set=vis)
        
        #Instead use for now MPRelaxSet for all by not specifying anything
        new_optimize_wf = OptimizeFW(structure)
        
        return FWAction(additions = [new_optimize_wf])
        
@explicit_serialize
class FT_PassSpec(FiretaskBase):
    """Update only certain keys in the first level of the spec.
    
    If the key_list contatins only '_all', update the whole spec!
    """
    
    _fw_name = 'Pass Spec'
    required_params=['key_list']
    def run_task(self, fw_spec):
        update = {}
        if self['key_list'] == ['_all']:
            spec = fw_spec
            return FWAction(update_spec = spec)
        else:
            for k in self['key_list']:
                if fw_spec.get(k) is None:
                    raise ValueError('{} can not be passed on to the next'
                                     'FT/FW because it is not in the spec.\n'
                                     'Currently the spec has the keys:\n'
                                     '{}'.format(k, fw_spec.keys()))
                update[k] = fw_spec.get(k)
            return FWAction(update_spec = update)

@explicit_serialize
class FT_CopyInSpec(FiretaskBase):
    """Copy something saved in the spec to another position in the spec.
    """
    
    _fw_name = 'Copy in Spec'
    required_params = ['in_loc', 'out_loc']
    def run_task(self, fw_spec):
        from HelperFunctions import GetValueFromNestedDict, SetInSpecFWAction
        thing = GetValueFromNestedDict(fw_spec, self['in_loc'])
        FWA = SetInSpecFWAction(self['out_loc'], thing)
        return FWA
        

@explicit_serialize
class FT_DoNothing(FiretaskBase):
    """Does Nothing
    """
    _fw_name = 'Do Nothing'
    def run_task(self, fw_spec):
        fw_spec['testlala'] ='test'

@explicit_serialize
class FT_PrintSpec(FiretaskBase):
    """Prints the spec of the current Workflow to the screen.
    
    Not only prints the current spec in a pretty way, but also returns a
    FWAction that updates the spec of future to include the current spec.
    """
    
    _fw_name = 'Print Spec'
    
    def run_task(self, fw_spec):
        import pprint
        
        pprint.pprint(fw_spec)
        
        spec = fw_spec
        return FWAction(update_spec = spec)        

@explicit_serialize
class FT_FetchStructureFromFormula(FiretaskBase):
    """Fetches a structure from the MP database and puts it in the spec.
    
    Uses the helper function GetLowEnergyStructure to get the structure for a
    given chemical formula with the lowest formation energy. If a mp_id is
    also given, this exact structure will be fetched. The structure will be 
    saved on the first level of the spec under the given name. Information on
    the metallicity of the material will also be saved in the spec.
    
    Parameters
    ----------
    materials_dict_loc: list of str
        Location where the formula of the structure to fetch is located.
        E.g. ['inputs', 'bottom_slab'] will lead to the formula located at
        fw_spec['inputs']['bottom']['formula'] to be loaded.
        The input key 'mp_id' will also be loaded and has to be present at
        the same location (e.g. fw_spec['inputs']['bottom']['mp-id']) but it
        may be None-type and then is not used to fetch the structure.
        
    structure_loc: list of str, optional
        If given saves the structure at this location in the fw_spec.
        E.g. ['my_structures', 'relaxed_structures', 'bulk_TiO'] will allow
        you to find the structure later in
        fw_spec[my_structures']['relaxed_structures']['bulk_TiO']
        If not given, the structure will be placed in the spec on the same
        level as the formula with the key: formula+'_fromMP'
    
    Returns
    -------
    An updated fw_spec that includes a new structure and updates the mp_id and
    information about the bandgap of the structure.
    """
    
    _fw_name = 'Fetch Structure From Spec'
    required_params = ['materials_dict_loc']
    optional_params = ['structure_loc']
    
    def run_task(self, fw_spec):
        from HelperFunctions import GetLowEnergyStructure, \
                                    GetValueFromNestedDict, GetGapFromMP, \
                                    UpdateNestedDict, WriteNestedDictFromList
        
        mat_dict = self['materials_dict_loc']
        formula_loc = mat_dict+['formula']
        mpid_loc = mat_dict+['mp_id']
        formula = GetValueFromNestedDict(fw_spec, formula_loc)
        mp_id = GetValueFromNestedDict(fw_spec, mpid_loc)
        
        
        structure, MP_ID = GetLowEnergyStructure(formula,
                                                 MP_ID=mp_id,
                                                 PrintInfo=True)
        
        bandgap = GetGapFromMP(MP_ID)
        if bandgap > 0.1:
            is_metal = False
        else:
            is_metal = True
            
        updated_mpid = WriteNestedDictFromList(mpid_loc, MP_ID, d={})
        updated_metal = WriteNestedDictFromList(mat_dict+['is_metal'],
                                                is_metal, d={})
        if 'structure_name' in self:
            updated_struct = WriteNestedDictFromList(self['structure_loc'],
                                                     structure, d={})
        else:
            updated_struct = WriteNestedDictFromList(mat_dict+[formula+
                                                               '_fromMP'],
                                                     structure, d={})
        
        
        import pprint
        spec = fw_spec
        spec = UpdateNestedDict(spec, updated_mpid)
        spec = UpdateNestedDict(spec, updated_metal)
        spec = UpdateNestedDict(spec, updated_struct)
        
        return FWAction(update_spec = spec)

@explicit_serialize
class FT_FetchStructure(FiretaskBase):
    """Fetches a structure from the MP database and puts it in the spec.
    
    Uses the helper function GetLowEnergyStructure to get the structure for a
    given chemical formula with the lowest formation energy. If a MP_ID is
    also given, this exact structure will be fetched. The 'structures' key of
    the spec is updated to include the fetched structur under the formula key.
    
    Parameters
    ----------
    formula: str
        Required input of a chemical formula
        e.g.: NaCl, Fe2O3, SiO, FeCW
    MP_ID: str
        Optional Input of Materials Project ID
        of the exact desired structure
        e.g. 'mp-990448'
    
    Yields
    ------
    updated structures dictionary in the spec.
    """
    
    _fw_name = 'Fetch Structure'
    required_params = ['formula', 'MP_ID']
    
    def run_task(self, fw_spec):
        from HelperFunctions import GetLowEnergyStructure
        
        structure, MP_ID = GetLowEnergyStructure(self['formula'],
                                                 MP_ID=self['MP_ID'],
                                                 PrintInfo=False)
        
        if 'structures' in fw_spec:
            structures = fw_spec['structures']
        else:
            structures = {}
            
        structures[self['formula']] = structure
        
        return FWAction(update_spec = structures)

        

@explicit_serialize
class FT_ReadInputFile(FiretaskBase):
    """Reads an input file and puts the contents in a dictionary.

    Reads data from an input file present in the execution directory and puts
    it in a dictionary. Empty lines and lines starting with '#' are ignored,
    as are trailing comments, which are designeated with a '#' character.
    The format of the input file has to follow the format:
        Flag_1 = Setting_1
        Flag_2 = Setting_2
        ...
        ...
        
    Parameters
    ----------
    filename: str
        File to be opened and read by the FireTask
    
    Returns
    -------
    dictionary is pushed to the spec and contains key-value pairs
        {'Flag_1': 'Setting_1',
         'Flag_2': 'Setting_2',
          ...}
    """
    
    _fw_name = 'Read input file'
    required_params = ['filename']
    
    def run_task(self, fw_spec):
        try:
            inputfile = open(self['filename'], 'r')
        except IOError:
            with open('output.txt', 'a') as out:
                out.write(' ')
                out.write('Error opening input file "'+self['filename']+
                          '". Please make sure that it is present.\n')
                out.write(' ')
            raise SystemExit
        #setup a dictionary to read the data into. This makes it possible to change the order of lines.
        input_dict={}
        for line in inputfile:
            #remove any leading blanks and check for comments or empty lines:
            line.strip()
            if not line.startswith('#') and not line.startswith('\n'):
                #remove comments and also trailing blanks
                line = line.partition('#')[0]
                line = line.rstrip()
                #partition the data in key and value pairs and put them in a dictionary
                name, var = line.partition('=')[::2]
                input_dict[name.strip()] = var.strip()

        inputfile.close()
        
        return FWAction(update_spec={self['filename']: input_dict})


@explicit_serialize
class FT_CheckCompParamDict(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values if needed.
    
    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.
    
    Parameters
    ----------
    input_dict_loc: list of str
        The location of a dictionary in the spec that holds the computational
        input parameters
    output_dict_loc: list of str, optional
        The location of the output dictionary that is going to be put into the
        fw_spec. If not specified this will default to the input_dict_loc.
    
    Returns
    -------
    dict
        A dictionary is pushed to the fw_spec. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys.
    """

    _fw_name = 'Check Input Dict'
    required_params = ['input_dict_loc']
    optional_params = ['output_dict_loc']
    
    def run_task(self, fw_spec):
        """Run the FireTask.

        Edit the essential and additional keys here if necessary, as well as
        the default values given for the additional keys.
        """
        from HelperFunctions import GetValueFromNestedDict, \
                                    WriteNestedDictFromList, \
                                    UpdateNestedDict
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['use_vdw']
        additional_keys = ['volume_tolerance', 'energy_tolerance', 'use_spin',
                           'BM_tolerance', 'functional']
        
        volume_tolerance_default = 0.001
        energy_tolerance_default = 0.001
        functional_default = 'SCAN'
        use_spin_default = True
        BM_tolerance_default = 0.01
        #####################################################################
        
        input_dict = GetValueFromNestedDict(fw_spec, self['input_dict_loc'])
        
        if 'output_dict_loc' in self:
            out_loc = self['output_dict_loc']
        else:
            out_loc = self['input_dict_loc']

        #Define all known keys here
        known_keys = essential_keys + additional_keys
        #Create lists for input options condidered to be True:
        true_list = ['true', 'True', 'TRUE', '.TRUE.', '.true.', True, 
                     'yes', 'Yes', 'YES', '.YES.', '.yes.']
        #initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
        for key in input_dict.keys():
        #Check for unknown parameters:
            if key not in known_keys:
                raise SystemExit('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n'
                              '{}'.format(known_keys))
            elif key == 'use_vdw':
                if input_dict[key] in true_list:
                    out_dict['use_vdw'] = True
                else:
                    out_dict['use_vdw'] = False
                check_essential[key] = True
        
        
        #check if all essential keys are present and print out missing ones.
        missing_keys=[]
        for key in check_essential.keys():
            if check_essential[key] == False:
                missing_keys.append(key)
        if not all(check_essential.values()):
            raise SystemExit('The essential input parameters "'+missing_keys+
                              '" are missing. Check your input file!\n')
        
        for key in additional_keys:
            if key == 'use_spin':
                if key in input_dict:
                    if input_dict[key] in true_list:
                        out_dict['use_spin'] = True
                    else:
                        out_dict['use_spin'] = False
                else:
                    out_dict['use_spin'] = use_spin_default
            if key == 'volume_tolerance':
                if key in input_dict:
                    out_dict['volume_tolerance'] = float(input_dict[key])
                else:
                    out_dict['volume_tolerance'] = volume_tolerance_default
            if key == 'energy_tolerance':
                if key in input_dict:
                    out_dict['energy_tolerance'] = float(input_dict[key])
                else:
                    out_dict['energy_tolerance'] = energy_tolerance_default
            if key == 'functional':
                if key in input_dict:
                    out_dict['functional'] = str(input_dict[key])
                else:
                    out_dict['functional'] = functional_default
            if key == 'BM_tolerance':
                if key in input_dict:
                    out_dict['BM_tolerance'] = float(input_dict[key])
                else:
                    out_dict['BM_tolerance'] = BM_tolerance_default
                
        spec = fw_spec
        final_output = WriteNestedDictFromList(out_loc, out_dict, d={})
        updated_spec = UpdateNestedDict(spec, final_output)
        return FWAction(update_spec=updated_spec)
    
@explicit_serialize
class FT_CheckInterfaceParamDict(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values if needed.
    
    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.
    
    Parameters
    ----------
    input_dict_loc: list of str
        The location of a dictionary in the spec that holds the computational
        input parameters
    output_dict_loc: list of str, optional
        The location of the output dictionary that is going to be put into the
        fw_spec. If not specified this will default to the input_dict_loc.
    
    Returns
    -------
    dict
        A dictionary is pushed to the fw_spec. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys.
    """

    _fw_name = 'Check Input Dict'
    required_params = ['input_dict_loc']
    optional_params = ['output_dict_loc']
    
    def run_task(self, fw_spec):
        """Run the FireTask.

        Edit the essential and additional keys here if necessary, as well as
        the default values given for the additional keys.
        """
        from HelperFunctions import GetValueFromNestedDict, \
                                    WriteNestedDictFromList, \
                                    UpdateNestedDict
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['max_area']
        additional_keys = ['interface_distance', 'max_mismatch',
                           'max_angle_diff', 'r1r2_tol']
        
        defaults = {'interface_distance': 2.0,
                    'max_mismatch': 0.01,
                    'max_angle_diff': 2.0,
                    'r1r2_tol': 0.05}
        #####################################################################
        
        input_dict = GetValueFromNestedDict(fw_spec, self['input_dict_loc'])
        
        if 'output_dict_loc' in self:
            out_loc = self['output_dict_loc']
        else:
            out_loc = self['input_dict_loc']

        #Define all known keys here
        known_keys = essential_keys + additional_keys
        #initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
        for key in input_dict.keys():
        #Check for unknown parameters:
            if key not in known_keys:
                raise SystemExit('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n'
                              .format(known_keys))
            elif key == 'max_area':
                out_dict['max_area'] = float(input_dict[key])
                check_essential[key] = True
        
        
        #check if all essential keys are present and print out missing ones.
        missing_keys=[]
        for key in check_essential.keys():
            if check_essential[key] == False:
                missing_keys.append(key)
        if not all(check_essential.values()):
            raise SystemExit('The essential input parameters "'+missing_keys+
                              '" are missing. Check your input file!\n')
        
        for key in additional_keys:
            if key in input_dict:
                out_dict[key] = float(input_dict[key])
            else:
                out_dict[key] = defaults[key]
                
        spec = fw_spec
        final_output = WriteNestedDictFromList(out_loc, out_dict, d={})
        updated_spec = UpdateNestedDict(spec, final_output)
        return FWAction(update_spec=updated_spec)

@explicit_serialize
class FT_CheckMaterialInputDict(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values if needed.
    
    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.
    
    Parameters
    ----------
    input_dict_name: list of str
        List of keys that specifiy the location of the input parameters
        dictionary in the spec
    output_dict_name: list of str, optional
        Location of the output dictionary that is going to be put into the
        spec. If not specified this will default to the input_dict_loc.
    
    Returns
    -------
    dict
        A dictionary is pushed to the spec of the workflow. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys. The name in the spec is output_dict_name.
    """

    _fw_name = 'Check Material Input Dict'
    required_params = ['input_dict_loc']
    optional_params = ['output_dict_loc']
    
    def run_task(self, fw_spec):
        """Run the FireTask."""
        from HelperFunctions import UpdateNestedDict, GetValueFromNestedDict, \
                                    WriteNestedDictFromList
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['formula', 'miller']
        additional_keys = ['min_vacuum', 'mp_id', 'min_thickness']
        
        MP_ID_default = None
        min_thickness_default = 10.0
        min_vacuum_default = 25.0
        #####################################################################
        input_dict = GetValueFromNestedDict(fw_spec, self['input_dict_loc'])
        if 'output_dict_name' in self:
            out_loc = self['output_dict_loc']
        else:
            out_loc = self['input_dict_loc']

        #Define all known keys here
        known_keys = essential_keys + additional_keys
        #initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
            
        for key in input_dict.keys():
        #Check for unknown parameters:
            if key not in known_keys:
                with open('output.txt', 'a') as out:
                    out.write(' ')
                    out.write('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n')
                    out.write(str(known_keys))
                    out.write('\n')
                raise SystemExit
            elif key == 'formula':
                out_dict['formula'] = str(input_dict[key])
                check_essential[key] = True
            elif key == 'miller':
                out_dict['miller'] = [int(k) for k in list(input_dict[key])]
                check_essential[key] = True
        
        #check if all essential keys are present and print out missing ones.
        for key in check_essential.keys():
            if check_essential[key] == False:
                with open('output.txt', 'a') as out:
                    out.write('The essential input parameter "'+key+
                              '" is missing. Check your input file!\n')
        if not all(check_essential.values()):
            with open('output.txt', 'a') as out:
                out.write('')
            raise SystemExit
            
        for key in additional_keys:
            if key == 'mp_id':
                if key in input_dict:
                    out_dict['mp_id'] = str(input_dict[key])
                else:
                    out_dict['mp_id'] = MP_ID_default
            if key == 'min_thickness':
                if key in input_dict:
                    out_dict['min_thickness'] = float(input_dict[key])
                else:
                    out_dict['min_thickness'] = min_thickness_default
            if key == 'min_vacuum':
                if key in input_dict:
                    out_dict['min_vacuum'] = float(input_dict[key])
                else:
                    out_dict['min_vacuum'] = min_vacuum_default
        
        spec = fw_spec
        spec_update = WriteNestedDictFromList(out_loc, out_dict, d={})
        updated_spec = UpdateNestedDict(spec, spec_update)
        return FWAction(update_spec=updated_spec)
            

@explicit_serialize
class FT_CheckInputDict(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values if needed.
    
    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.
    
    Parameters
    ----------
    input_dict_name: str
        The name of a dictionary that was previously constructed from an input
        file of the same name.
    output_dict_name: str, optional
        Name of the output dictionary that is going to be put into the workflow
        spec. If not specified this will default to the input_dict_name.
    
    Returns
    -------
    dict
        A dictionary is pushed to the spec of the workflow. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys. The name in the spec is output_dict_name.
    """

    _fw_name = 'Check Input Dict'
    required_params = ['input_dict_name']
    optional_params = ['output_dict_name']
    
    def run_task(self, fw_spec):
        """Run the FireTask.

        Edit the essential and additional keys here if necessary, as well as
        the default values given for the additional keys.
        """
        
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['formula', 'miller', 'use_vdw', 'use_spin']
        additional_keys = ['volume_tolerance', 'energy_tolerance',
                             'functional', 'BM_tolerance', 'MP_ID',
                             'min_thickness']
        
        volume_tolerance_default = 0.001
        energy_tolerance_default = 0.001
        functional_default = 'SCAN'
        BM_tolerance_default = 0.01
        MP_ID_default = None
        min_thickness_default = 10.0
        #####################################################################
        
        input_dict = fw_spec[self['input_dict_name']]
        if 'output_dict_name' in self:
            out_name = self['output_dict_name']
        else:
            out_name = self['input_dict_name']

        #Define all known keys here
        known_keys = essential_keys + additional_keys
        #initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
        for key in input_dict.keys():
        #Check for unknown parameters:
            if key not in known_keys:
                with open('output.txt', 'a') as out:
                    out.write(' ')
                    out.write('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n')
                    out.write(str(known_keys))
                    out.write('\n')
                raise SystemExit
            elif key == 'formula':
                out_dict['formula'] = str(input_dict[key])
                check_essential[key] = True
            elif key == 'miller':
                out_dict['miller'] = [int(k) for k in list(input_dict[key])]
                check_essential[key] = True
            elif key == 'use_vdw':
                out_dict['use_vdW'] = bool(input_dict[key])
                check_essential[key] = True
            elif key == 'use_spin':
                out_dict['use_spin'] = bool(input_dict[key])
                check_essential[key] = True
        
        
        #check if all essential keys are present and print out missing ones.
        for key in check_essential.keys():
            if check_essential[key] == False:
                with open('output.txt', 'a') as out:
                    out.write('The essential input parameter "'+key+
                              '" is missing. Check your input file!\n')
        if not all(check_essential.values()):
            with open('output.txt', 'a') as out:
                out.write('')
            raise SystemExit
        
        for key in additional_keys:
            if key == 'volume_tolerance':
                if key in input_dict:
                    out_dict['volume_tolerance'] = float(input_dict[key])
                else:
                    out_dict['volume_tolerance'] = volume_tolerance_default
            if key == 'energy_tolerance':
                if key in input_dict:
                    out_dict['energy_tolerance'] = float(input_dict[key])
                else:
                    out_dict['energy_tolerance'] = energy_tolerance_default
            if key == 'functional':
                if key in input_dict:
                    out_dict['functional'] = str(input_dict[key])
                else:
                    out_dict['functional'] = functional_default
            if key == 'MP_ID':
                if key in input_dict:
                    out_dict['MP_ID'] = str(input_dict[key])
                else:
                    out_dict['MP_ID'] = MP_ID_default
            if key == 'BM_tolerance':
                if key in input_dict:
                    out_dict['BM_tolerance'] = float(input_dict[key])
                else:
                    out_dict['BM_tolerance'] = BM_tolerance_default
            if key == 'min_thickness':
                if key in input_dict:
                    out_dict['min_thickness'] = float(input_dict[key])
                else:
                    out_dict['min_thickness'] = min_thickness_default
                
        return FWAction(update_spec={out_name: out_dict})
        
@explicit_serialize
class FT_MakeHeteroStructure(FiretaskBase):
    """Matches two slab systems to form a heterogeneous interface.

    If the match fails, the accuracy criteria are relaxed in steps of 5% until
    a match is found.

    Parameters
    ----------
    bottom_slab_loc: list of str
        list of keys in the spec that point to the bottom slab.
    top_slab_loc: list of str
        list of keys in the spec that point to the top slab.
    parameters_loc: list of str
        location of the dictionary in the spec holding the following keys:
        interface_distance (float): Distance between the two materials.
        max_area (float): Maximal cross section area of the matched cell in
            Angstrom squared. Defaults to 200.,
        max_mismatch (float): Maximal allowed mismatch between lattice vector
            length. Defaults to 0.01 (1%),
        max_angle_diff (float): Maximal allowed mismatch between lattice
            vector angle in degrees. Defaults to 1.
        r1r2_tol (float): Tolerance between factors r1 and r2 which relate to
            the matched cell areas as:
            abs(r1*area1-r2*area2)/max_area <= r1r2_tol. Defaults to 0.02
    out_loc: list of str, optional
        If given saves the output structure, and the two matched slabs in the
        fw_spec at the location specified. If not given the location defaults
        to: ['interface'], where the slabs and interfaces will be saved as:
            'top_slab'
            'bottom_slab'
            'initial_match'

    Returns
    -------
    Matched interface structure and bottom and top slabs in the fw_spec.
    """
    
    _fw_name = 'Make Hetero Structure'
    required_params = ['bottom_slab_loc', 'top_slab_loc', 'parameters_loc']
    optional_params = ['out_loc']
    
    def run_task(self, fw_spec):
        from mpinterfaces.transformations import get_aligned_lattices
        from mpinterfaces.transformations import generate_all_configs
        from HelperFunctions import UpdateNestedDict, \
                                    GetValueFromNestedDict, \
                                    WriteNestedDictFromList
        parameters = GetValueFromNestedDict(fw_spec, self['parameters_loc'])
        bottom_slab = GetValueFromNestedDict(fw_spec, self['bottom_slab_loc'])
        top_slab = GetValueFromNestedDict(fw_spec, self['top_slab_loc'])
        
        bottom_aligned, top_aligned = get_aligned_lattices(
                bottom_slab,
                top_slab,
                max_area = parameters['max_area'],
                max_mismatch = parameters['max_mismatch'],
                max_angle_diff = parameters['max_angle_diff'],
                r1r2_tol = parameters['r1r2_tol'])
        
        if bottom_aligned is not None:
# ============================================================================
# TODO: Find out if this is actually useful for the PES to return
#       all the interface structures we need, or if another stacking function
#       should be written with another input of lateral shifts.
#       Check out AdsorbateSiteFinder of pymatgen.analysis.adsorption!!
#       If generate all configs will work, we just have to remove the
#       [0] index and move the whole list of structures to the spec.
# ============================================================================
            hetero_interfaces = generate_all_configs(top_aligned,
                                                     bottom_aligned,
                                                     nlayers_2d=1,
                                                     nlayers_substrate=1,
                            seperation=parameters['interface_distance'])
        
            name = self['bottom_slab_loc'][-1]+self['bottom_slab_loc'][-1]
            hetero_interfaces[0].to(fmt='poscar', filename=
                                    'POSCAR_'+name+'.vasp')
            
            if 'out_loc' in self:
                out_loc = self['out_loc']
            else:
                out_loc = ['interface']
            
            out_dict = WriteNestedDictFromList(out_loc,
                                        {'top_slab': top_aligned,
                                         'bottom_slab': bottom_aligned,
                                         'intial_match': hetero_interfaces[0]},
                                        d={})
            
            spec = fw_spec
            updated_spec = UpdateNestedDict(spec, out_dict)
            return FWAction(update_spec=updated_spec)
        else:
            f = 1.05 # factor with wich to multiply the missmatch criteria
            new_parameters = {'max_area': parameters['max_area'] * f,
                    'max_mismatch': parameters['max_mismatch'] * f,
                    'max_angle_diff': parameters['max_angle_diff'] * f,
                    'r1r2_tol': parameters['r1r2_tol'] * f,
                    'interface_distance': parameters['interface_distance']
                                }
            spec = fw_spec
            parameters_dict = WriteNestedDictFromList(self['parameters_loc'],
                                                      new_parameters, d={})
            updated_spec = UpdateNestedDict(spec, parameters_dict)
            new_fw = Firework(FT_MakeHeteroStructure(
                            bottom_slab_name=self['bottom_slab_name'],
                            top_slab_name=self['top_slab_name'],
                            parameters_loc=self['parameters_loc']),
                            spec=updated_spec,
                            name = 'Make Heterostructure with relaxed params')
            return FWAction(detours = new_fw)

@explicit_serialize
class FT_MakeSlabFromStructure(FiretaskBase):
    """Makes a slab with certain orientation out of a bulk structure.

    The slab has a defined thickness and orientation and the bulk structure is
    loaded from the fw_spec

    Parameters
    ----------
    bulk_loc: list of str
        list of keys in the spec that point to the bulk input structure.
    dict_loc: list of str
        list of keys for the dictionary in the spec that contains input data
        in the following keys:
            miller: list of int
            min_thickness: float
            min_vacuum: float
    out_loc: list of str, optional
        If given saves the output structure at this location in the fw_spec.
        E.g. ['my_structures', 'unrelaxed_structures', 'TiO111_slab'] will
        allow you to find the structure later in
        fw_spec[my_structures']['unrelaxed_structures']['TiO111_slab']
        If not given, the structure will be placed in the spec on the same
        level as the input structure with the key: formula+miller
    
    Returns
    -------
        Pymatgen slab structure to the spec
    """
    
    _fw_name = 'Make slab from structure in spec'
    required_params = ['bulk_loc', 'dict_loc']
    optional_params = ['out_loc']
    
    def run_task(self, fw_spec):
        from mpinterfaces.interface import Interface
        from HelperFunctions import UpdateNestedDict, GetValueFromNestedDict,\
                                    WriteNestedDictFromList
        
        bulk_structure = GetValueFromNestedDict(fw_spec, self['bulk_loc'])
        param_dict = GetValueFromNestedDict(fw_spec, self['dict_loc'])
            
        miller = param_dict['miller']
        min_thickness = param_dict['min_thickness']
        min_vacuum = param_dict['min_vacuum']

        
        #Construct the slab out of the bulk_structure
        slab = Interface(bulk_structure, hkl = miller,
                   min_thick = min_thickness,
                   min_vac = min_vacuum, primitive = False,
                   from_ase = True)

        spec = fw_spec
        miller = ''.join(str(e) for e in miller)
        slab_name = bulk_structure.composition.reduced_formula + miller
        if 'out_loc' in self:
            out_loc = self['out_loc']
        else:
            out_loc = self['bulk_loc'][:-1]+[slab_name]
        struct_dict = WriteNestedDictFromList(out_loc, slab, d={})
        updated_spec = UpdateNestedDict(spec, struct_dict)
        return FWAction(update_spec=updated_spec)

@explicit_serialize
class FT_MakeSlabFromFormula(FiretaskBase):
    """Makes a slab with certain orientation out of a bulk structure.

    The slab has a defined thickness and the bulk structure is loaded from the
    MaterialsProject database.

    Args
    ----
        formula (str):      Chemical formula of the structure (e.g. 'Fe2O3').
                            The structure with the lowest formation energy for
                            this formula is used unless MP_ID is specified.
        mp_id (str):        Materials Project ID of the selected bulk structure
                            (e.g. 'mp-990448').
        miller (list):      List of miller indices to describe the selected
                            surface of the slab (e.g. [1, 1, 0]).
        thickness (float):  Minimal thickness for the slab in Angstrom.
        vacuum (float):     Thickness of the vacuum layer between repeating
                            slabs. Defaults to 25.0 Angstrom
    ----
    
    Returns
    -------
        Pymatgen slab structure
    ------- 
    """
    
    _fw_name = 'Make Slab from Formula'
    required_params = ['formula', 'miller', 'thickness']
    optional_params = ['mp_id', 'vacuum']
    
    def run_task(self, fw_spec):
        from HelperFunctions import GetLowEnergyStructure
        from mpinterfaces.interface import Interface
        
        if 'mp_id' in self:
            MP_ID = self['mp_id']
        else:
            MP_ID = None
        
        if 'vacuum' in self:
            vacuum = self['vacuum']
        else:
            vacuum = 25.0
            
        #Get bulk structure from MP database
        bulk_structure, MP_ID_bulk = GetLowEnergyStructure(self['formula'],
                                               MP_ID=MP_ID,
                                               PrintInfo=False)

        
        #Construct the slab out of the bulk_structure
        slab = Interface(bulk_structure, hkl = self['miller'],
                   min_thick = self['thickness'],
                   min_vac = vacuum, primitive = False,
                   from_ase = True)
        
        miller = ''.join(str(e) for e in self['miller'])
        slab_name = self['formula'] + miller
        return FWAction(update_spec={slab_name: slab}, 
                        stored_data={slab_name: slab})

@explicit_serialize
class FT_WritePrecalc(FiretaskBase):
    """Write a PRECALC file for the getKPoints script."""
    
    _fw_name = "Write PRECALC file"
    required_params = ['PrecalcDict']
    def run_task(self, fw_spec):
        from HelperFunctions import WriteFileFromDict
        WriteFileFromDict(self['PrecalcDict'], 'PRECALC')

@explicit_serialize
class FT_WriteGeneralKpointsInputs(FiretaskBase):
    """Prepares the VASP-4 inputs for the generation of a generalized K-mesh.

    Prepares the necessary files (POSCAR, PRECALC, maybe INCAR)
    for the K-Point Grid Generator of the Mueller group at John Hopkins
    http://muellergroup.jhu.edu/K-Points.html

    Args
    ----
        structure (pymatgen structure object): Required input structure for
                                               which the Kpoints mesh should
                                               be created. 
        PrecalcDict (dict): Required parameters for the Grid Generator, should
                            definately include the MINDISTANCE flag.
        IncarDict (dict):   Optional INCAR file as a dictionary to pass
                            information about ISYM, MAGMOM, etc...
        WorkDir (str):      Optional working directory where everything happens
    ----
    """
    
    _fw_name = "Write Generalized Kpoints Inputs"
    required_params = ['structure', 'PrecalcDict']
    optional_params = ['IncarDict', 'WorkDir']

    def run_task(self, fw_spec):
        from HelperFunctions import WriteFileFromDict
        
        if 'WorkDir' in self:
            workdir = self['WorkDir']
        else:
            workdir = './'
        
        self['structure'].to(fmt='poscar', filename=workdir+'POSCAR')
        WriteFileFromDict(self['PrecalcDict'], workdir+'PRECALC')

        if 'IncarDict' in self:
            WriteFileFromDict(self['IncarDict'], workdir+'INCAR')
            

@explicit_serialize
class FT_GetKpoints(FiretaskBase):
    """Reads a KPOINTS file and pushes it to the spec.

    Args
    ----
        workdir (str): Optional working directory. Defaults
                       to './'
    ----
    
    Returns
    -------
        FWAction
    -------
    """
    
    _fw_name = 'Get Kpoints'
    optional_params = ['WorkDir', 'KPTS_loc']

    def run_task(self, fw_spec):
        import os
        from pymatgen.io.vasp.inputs import Kpoints
        from HelperFunctions import WriteNestedDictFromList
        
        if 'WorkDir' in self:
            workdir = self['WorkDir']
        else:
            workdir = './'
        

        #Make a pymatgen Kpoints Object from a KPOINTS file
        KPTS = Kpoints().from_file(workdir+'KPOINTS')
        
        if 'KPTS_loc' in self:
            out_dict = WriteNestedDictFromList(self['KPTS_loc'], KPTS, d={})
        else:
            out_dict = {'KPOINTS': KPTS}

        #Check for 'KPOINTS.log in workdir and read if it is there'
        if os.path.isfile(workdir+'KPOINTS.log'):
            with open(workdir+'KPOINTS.log') as f:
                Kpoints_Info = f.readlines()
            return FWAction(stored_data={'KPOINTS_Info': Kpoints_Info},
                             mod_spec=[{'_push': out_dict}])
        else:
            return FWAction(mod_spec=[{'_push': out_dict}])


@explicit_serialize
class FT_CleanUpGeneralizedKpointFiles(FiretaskBase):
    """Cleans up files that are no longer needed after K-mesh generation.
    
    Cleans up all the files needed and created by the 
    K-Point Grid Generator of the Mueller group at John Hopkins
    http://muellergroup.jhu.edu/K-Points.html
    """
    
    _fw_name = 'Clean Up Generalized Kpoint Files'
    optional_params = ['WorkDir']

    def run_task(self, fw_spec):
        from HelperFunctions import RemoveMatchingFiles
        
        if self['WorkDir']:
            workdir = self['WorkDir']
        else:
            workdir = './'

        RemoveMatchingFiles([workdir+'KPOINTS*',
                            workdir+'POSCAR*',
                            workdir+'INCAR', 
                            workdir+'PRECALC'])

