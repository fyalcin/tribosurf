#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:30:31 2020

@author: mwo
"""
from fireworks import FWAction, FiretaskBase, Firework, Workflow, LaunchPad
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.utils.utils import load_class
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.config import HALF_KPOINTS_FIRST_RELAX, RELAX_MAX_FORCE, \
    VASP_CMD, DB_FILE
    
from HelperFunctions import *

@explicit_serialize
class WriteVaspFromIOSet_Mod(FiretaskBase):
    """
    Create VASP input files using implementations of pymatgen's AbstractVaspInputSet. An input set
    can be provided as an object or as a String/parameter combo.

    Required params:
        structure_spec_loc (list of str): Location of the input structure
                in the fw_spec. E.g. ['structures', 'to_relax', 'MyStruct']
                will point to fw_spec['structures']['to_relax']['MyStruct']
        vasp_input_set (AbstractVaspInputSet or str): Either a VaspInputSet object or a string
            name for the VASP input set (e.g., "MPRelaxSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set, use this as a dict
            to specify kwargs for instantiating the input set parameters. For example, if you want
            to change the user_incar_settings, you should provide: {"user_incar_settings": ...}.
            This setting is ignored if you provide the full object representation of a VaspInputSet
            rather than a String.
    """

    required_params = ["structure_spec_loc", "vasp_input_set"]
    optional_params = ["vasp_input_params"]

    def run_task(self, fw_spec):
        structure = GetValueFromNestedDict(fw_spec, self['structure_spec_loc'])
        # if a full VaspInputSet object was provided
        if hasattr(self['vasp_input_set'], 'write_input'):
            vis = self['vasp_input_set']

        # if VaspInputSet String + parameters was provided
        else:
            vis_cls = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])
            vis = vis_cls(structure, **self.get("vasp_input_params", {}))
        vis.write_input(".")



class OptimizeFW_Mod(Firework):

    def __init__(self, structure_spec_loc,
                 name="structure optimization", vasp_input_set='MPRelaxSet',
                 vasp_cmd=VASP_CMD, override_default_vasp_params=None,
                 ediffg=None, db_file=DB_FILE,
                 force_gamma=True, job_type="double_relaxation_run",
                 max_force_threshold=RELAX_MAX_FORCE,
                 auto_npar=">>auto_npar<<",
                 half_kpts_first_relax=HALF_KPOINTS_FIRST_RELAX, parents=None,
                 **kwargs):
        """
        Optimize the given structure.

        Args:
            structure_spec_loc (list of str): Location of the input structure
                in the fw_spec. E.g. ['structures', 'to_relax', 'MyStruct']
                will point to fw_spec['structures']['to_relax']['MyStruct']
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            ediffg (float): Shortcut to set ediffg in certain jobs
            db_file (str): Path to file specifying db credentials to place output parsing.
            force_gamma (bool): Force gamma centered kpoint generation
            job_type (str): custodian job type (default "double_relaxation_run")
            max_force_threshold (float): max force on a site allowed at end; otherwise, reject job
            auto_npar (bool or str): whether to set auto_npar. defaults to env_chk: ">>auto_npar<<"
            half_kpts_first_relax (bool): whether to use half the kpoints for the first relaxation
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        
        override_default_vasp_params = override_default_vasp_params or {}

        t = []
        t.append(WriteVaspFromIOSet_Mod(structure_spec_loc=structure_spec_loc,
                                    vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type,
                                  max_force_threshold=max_force_threshold,
                                  ediffg=ediffg,
                                  auto_npar=auto_npar,
                                  half_kpts_first_relax=half_kpts_first_relax))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(OptimizeFW_Mod, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             structure_spec_loc[-1], name),
                                         **kwargs)
    

def GetLowEnergyStructure(chem_formula, MP_ID=None, PrintInfo=False):
    """
    A function that searches the MaterialsProject Database
    for structures that match the given chemical formula
    and selcts the one with the lowest formation energy
    per atom. If several 
    Inputs: 
    chem_formula (str): Required input of a chemical formula
                        e.g.: NaCl, Fe2O3, SiO, FeCW
    MP_ID (str):        Optional Input of Materials Project ID
                        of the exact desired structure
                        e.g. 'mp-990448'
    PrintInfo (bool):   Optional variable defaults to "False".
                        if set to "True", will print some
                        information about how many structures
                        have been found and the ID of the selcted
                        one.
                        
    Returns: pymatgen Structure object and associated MP_ID
    """
    from pymatgen import MPRester

    # load structure from MaterialsProjct Website
    with MPRester() as mpr:
        if MP_ID:
            struct = mpr.get_structure_by_material_id(MP_ID)
            return struct, MP_ID
        else:
            id_list = mpr.query(criteria={'pretty_formula': chem_formula,
                                        'e_above_hull': 0.0},
                              properties=['material_id'])
            if id_list == []:
                raise NameError('{} has not been found in the MaterialsProject'
                            'database'.format(chem_formula))
            else:
                MP_ID = id_list[0]['material_id']
                struct = mpr.get_structure_by_material_id(MP_ID)
                return struct, MP_ID

def GetValueFromNestedDict(dictionary,key_list):
    """Take a list of keys and a nested dictionary and return the value.
    
    A (usually nested) dictionary is given along side a list of keys to this
    dictionary. The value at the end of this key_list is returned using a
    recursive method. E.g. key_list=['key_a', 'key_b', 'key_c'] will return
    dictionary['key_a']['key_b']['key_c']
    

    Parameters
    ----------
    dictionary : dict
        Input dictionary which can be nested
    key_list : list of str
        The list of keys to the nested dictionary at the end of which the
        value that is needed is located

    Returns
    -------
        Returns whatever is stored at the appropriate position in the
        dictionary. Type can vary!

    """
    if len(key_list) > 1:
        return GetValueFromNestedDict(dictionary[key_list[0]], key_list[1:])
    return dictionary[key_list[0]]


@explicit_serialize
class FT_FetchStructureFromInput(FiretaskBase):
    """Fetches a structure from the MP database and puts it in the spec.
    
    Uses the helper function GetLowEnergyStructure to get the structure for a
    given chemical formula with the lowest formation energy. If a MP_ID is
    also given, this exact structure will be fetched. The 'structures' key of
    the spec is updated to include the fetched structur under the formula key.
    
    Parameters
    ----------
    input_dict_name: str
        Required input dictionary where a chemical formula for the structure
        is given under the key 'formula'. The input key 'MP_ID' must also  be
        present but might be 'None' and will then not be used.
    ---------
    
    Yields
    ------
    updated structures dictionary in the spec.
    ------
    """
    
    _fw_name = 'Fetch Structure From Input'
    required_params = ['input_dict_name']
    
    def run_task(self, fw_spec):
        input_dict = fw_spec[self['input_dict_name']]
        
        if 'mp_id' in input_dict:
            mp_id = input_dict['mp_id']
        else:
            mp_id = None
        
        structure, MP_ID = GetLowEnergyStructure(input_dict['formula'],
                                                 MP_ID=mp_id,
                                                 PrintInfo=False)
        
        if 'structures' in fw_spec:
            structures = fw_spec['structures']
        else:
            structures = {}          
        structures[input_dict['formula']] = structure
        
        spec = fw_spec
        spec.update({'structures': structures})
        return FWAction(update_spec = spec)



material = {'formula': 'Fe', 'miller': [0, 0, 1], 'min_thickness': 10}

fw_1 = Firework(FT_FetchStructureFromInput(input_dict_name='material'),
                spec={'material': material})

fw_2 = OptimizeFW_Mod(structure_spec_loc=['structures', material['formula']],
                      parents=[fw_1])
wf = Workflow([fw_1, fw_2])

# finally, instatiate the LaunchPad and add the workflow to it
lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(wf)

rapidfire(lpad)