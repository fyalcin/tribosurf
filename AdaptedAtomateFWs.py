# coding: utf-8

import warnings



"""
Standard Fireworks from the atomate package that have been adapted to read
structures from the spec and write output to the spec.
"""

from fireworks import Firework
from fireworks import FiretaskBase, explicit_serialize

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
        t.append(VaspToDb(db_file=db_file,
                          additional_fields={"task_label": name}))
        super(OptimizeFW_Mod, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             structure_spec_loc[-1], name),
                                         **kwargs)
