""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from pymatgen.io.vasp.inputs import Poscar
from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

@explicit_serialize
class FT_MakeHeteroStructure(FiretaskBase):
    """Matches two slab systems to form a heterogeneous interface.

    If the match fails, the accuracy criteria are relaxed in steps of 5% until
    a match is found.

    Returns
    -------
    Matched interface structure and bottom and top slabs in the fw_spec.
    """
    
    _fw_name = 'Make Hetero Structure'
    required_params = ['mp_id_1', 'miller_1', 'mp_id_2', 'miller_2']
    optional_params = ['db_file']
    
    def run_task(self, fw_spec):

        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        miller_1 = self.get('miller_1')
        miller_2 = self.get('miller_2')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
            
        db = GetHighLevelDB(db_file)
        slab_col = 
        inter_col = db[functional+'.interface_data']
        
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



