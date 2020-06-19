""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from pymatgen.io.vasp.inputs import Poscar
from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

   
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



