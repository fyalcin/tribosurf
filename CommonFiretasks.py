#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python
# -*- coding: utf-8 -*-
"""A collection of Firetasks to be used for the FireFlow project.

Created on Mon Mar 23 14:33:47 2020
@author: mwo
"""
from fireworks.core.firework import FWAction, FiretaskBase, Firework
from fireworks.utilities.fw_utilities import explicit_serialize

# =============================================================================
# Custom FireTasks
# =============================================================================


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
class FT_PrintSpec(FiretaskBase):
    """Prints the spec of the current Workflow to the screen 
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
    saved on the first level of the spec under the given name.
    
    Parameters
    ----------
    formula_loc: list of str
        Location where the formula of the structure to fetch is located.
        E.g. ['inputs', 'bottom_slab', 'formula']
        will point to fw_spec['inputs']['bottom']['formula']
        If there is an input key 'mp_id' at the same level as the formula
        (e.g. fw_spec['inputs']['bottom']['mp-id']) it will be used to
        exactly define the structure to be fetched.
        
    structure_name: str
        If given saves the structrue under this name in the first level of the
        fw_spec. E.g. fw_spec[structure_name] = structure
    ---------
    
    Returns
    ------
    An updated fw_spec that includes a new structure.
    ------
    """
    
    _fw_name = 'Fetch Structure From Spec'
    required_params = ['formula_loc', 'structure_name']
    
    def run_task(self, fw_spec):
        from HelperFunctions import GetLowEnergyStructure, \
                                    GetValueFromNestedDict
        
        formula = GetValueFromNestedDict(fw_spec, self['formula_loc'])
        mp_id = GetValueFromNestedDict(fw_spec, 
                                       self['formula_loc'][:-1]+['mp_id'])
        
        structure, MP_ID = GetLowEnergyStructure(formula,
                                                 MP_ID=mp_id,
                                                 PrintInfo=False)
        
        spec = fw_spec
        spec.update({self['structure_name']: structure})
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
    ---------
    
    Yields
    ------
    updated structures dictionary in the spec.
    ------
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
        
    Args
    ----
        filename (str):     File to be opened and read by the FireTask
    ----
    
    Yields
    ------
        dictionary:         Is pushed to the spec and contains key-value
                            pairs {'Flag_1': 'Setting_1',
                                   'Flag_2': 'Setting_2',
                                   ...} 
    ------
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
class FT_CheckMaterialInputDict(FiretaskBase):
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
    ----------
    
    Returns
    -------
    dict
        A dictionary is pushed to the spec of the workflow. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys. The name in the spec is output_dict_name.
    -------
    """

    _fw_name = 'Check Material Input Dict'
    required_params = ['input_dict_name']
    optional_params = ['output_dict_name']
    
    def run_task(self, fw_spec):
        """Run the FireTask."""
        from HelperFunctions import UpdateNestedDict
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['formula', 'miller']
        additional_keys = ['min_vacuum', 'mp_id', 'min_thickness']
        
        MP_ID_default = None
        min_thickness_default = 10.0
        min_vacuum_default = 25.0
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
        updated_spec = UpdateNestedDict(spec, {out_name: out_dict})
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
    ----------
    
    Returns
    -------
    dict
        A dictionary is pushed to the spec of the workflow. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys. The name in the spec is output_dict_name.
    -------
    """

    _fw_name = 'Check Input Dict'
    required_params = ['input_dict_name']
    optional_params = ['output_dict_name']
    
    def run_task(self, fw_spec):
        """Run the FireTask

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

    Args
    ----
        bottom_slab_name (str): Name of the structure of the bottom slab that
                                needs to be aligned (structure itself must be
                                in fw_spec).
        top_slab_name (str):    Name of the structure of the top slab that
                                needs to be aligned (structure itself must be
                                in fw_spec).
        parameters (dict):  Dictionary containing the key (values):
                            interface_distance (float): Distance between the
                                two materials.
                            max_area (float): Maximal cross section
                                area of the matched cell in Angstrom squared.
                                Defaults to 200.,
                            max_mismatch (float): Maximal allowed mismatch
                                between lattice vector length in %.
                                Defaults to 0.01,
                            max_angle_diff (float): Maximal allowed mismatch
                                between lattice vector angle in degrees.
                                Defaults to 1.
                            r1r2_tol (float): Tolerance between factors r1 and
                                r2 which relate to the matched cell areas as:
                                abs(r1*area1-r2*area2)/max_area <= r1r2_tol.
                                Defaults to 0.02
    ----

    Returns
    -------
        Matched Structure (pymatgen structure object)
    -------
        
    """
    
    _fw_name = 'Make Hetero Structure'
    required_params = ['bottom_slab_name', 'top_slab_name', 'parameters']
    
    def run_task(self, fw_spec):
        from mpinterfaces.transformations import get_aligned_lattices
        from mpinterfaces.transformations import generate_all_configs
        
        bottom_aligned, top_aligned = get_aligned_lattices(
                fw_spec[self['bottom_slab_name']],
                fw_spec[self['top_slab_name']],
                max_area = self['parameters']['max_area'],
                max_mismatch = self['parameters']['max_mismatch'],
                max_angle_diff = self['parameters']['max_angle_diff'],
                r1r2_tol = self['parameters']['r1r2_tol'])
        
        if bottom_aligned is not None:
            #TODO: Find out if this is actually useful for the PES to return
            #all the interface structures we need, or if another stacking function
            #should be written with another input of lateral shifts.
            hetero_interfaces = generate_all_configs(top_aligned, bottom_aligned,
                            nlayers_2d=1, nlayers_substrate=1,
                            seperation=self['parameters']['interface_distance'])
        
            hetero_name = self['bottom_slab_name']+self['top_slab_name']
            hetero_interfaces[0].to(fmt='poscar', filename=
                                    'POSCAR_'+hetero_name+'.vasp')
        
            return FWAction(update_spec={'matched_interface':
                                         hetero_interfaces[0]}, 
                            stored_data={'matched_interface':
                                         hetero_interfaces[0]})
        else:
            f = 1.05 # factor with wich to multiply the missmatch criteria
            new_parameters = {'max_area': self['parameters']['max_area'] * f,
                    'max_mismatch': self['parameters']['max_mismatch'] * f,
                    'max_angle_diff': self['parameters']['max_angle_diff'] * f,
                    'r1r2_tol': self['parameters']['r1r2_tol'] * f,
                    'interface_distance': self['parameters']['interface_distance']
                                }
            new_fw = Firework(FT_MakeHeteroStructure(
                            bottom_slab_name=self['bottom_slab_name'],
                            top_slab_name=self['top_slab_name'],
                            parameters=new_parameters), fw_spec)
            return FWAction(additions = new_fw)

@explicit_serialize
class FT_MakeSlabFromStructure(FiretaskBase):
    """Makes a slab with certain orientation out of a bulk structure.

    The slab has a defined thickness and orientation and the bulk structure is
    loaded from the fw_spec

    Parameters
    ----------
    bulk_name: str
        name of the bulk structure in the spec from which the slab is made
    dict_name: str
        name of the dictionary in the spec that contains input data in the
        following keys:
            miller: list of int
            min_thickness: float
            min_vacuum: float
    ----
    
    Returns
    -------
        Pymatgen slab structure to the spec
    ------- 
    """
    
    _fw_name = 'Make slab from structure in spec'
    required_params = ['bulk_name', 'dict_name']
    
    def run_task(self, fw_spec):
        from mpinterfaces.interface import Interface
        
        bulk_structure = fw_spec[self['bulk_name']]
            
        miller = fw_spec[self['dict_name']]['miller']
        min_thickness = fw_spec[self['dict_name']]['min_thickness']
        min_vacuum = fw_spec[self['dict_name']]['min_vacuum']

        
        #Construct the slab out of the bulk_structure
        slab = Interface(bulk_structure, hkl = miller,
                   min_thick = min_thickness,
                   min_vac = min_vacuum, primitive = False,
                   from_ase = True)
        
        miller = ''.join(str(e) for e in miller)
        slab_name = bulk_structure.composition.reduced_formula + miller
        print(slab_name)
        print(slab_name)
        print(slab_name)
        return FWAction(update_spec={slab_name: slab}, 
                        stored_data={slab_name: slab})

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
        MP_ID (str):        Materials Project ID of the selected bulk structure
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
    optional_params = ['MP_ID', 'vacuum']
    
    def run_task(self, fw_spec):
        from HelperFunctions import GetLowEnergyStructure
        from mpinterfaces.interface import Interface
        
        if self['MP_ID']:
            MP_ID = self['MP_ID']
        else:
            MP_ID = None
        
        if self['vacuum']:
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
        
        if self['WorkDir']:
            workdir = self['WorkDir']
        else:
            workdir = './'
        
        self['structure'].to(fmt='poscar', filename=workdir+'POSCAR')
        WriteFileFromDict(self['PrecalcDict'], workdir+'PRECALC')

        if self['IncarDict']:
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
    optional_params = ['WorkDir']

    def run_task(self, fw_spec):
        import os
        from pymatgen.io.vasp.inputs import Kpoints
        
        if self['WorkDir']:
            workdir = self['WorkDir']
        else:
            workdir = './'

        #Make a pymatgen Kpoints Object from a KPOINTS file
        KPTS = Kpoints().from_file(workdir+'KPOINTS')

        #Check for 'KPOINTS.log in workdir and read if it is there'
        if os.path.isfile(workdir+'KPOINTS.log'):
            with open(workdir+'KPOINTS.log') as f:
                Kpoints_Info = f.readlines()
            return FWAction(stored_data={'KPOINTS_Info': Kpoints_Info},
                             mod_spec=[{'_push': {'KPOINTS': KPTS}}])
        else:
            return FWAction(mod_spec=[{'_push': {'KPOINTS': KPTS}}])


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

