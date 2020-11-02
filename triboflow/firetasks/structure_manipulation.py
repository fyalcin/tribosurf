""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""
import monty
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from fireworks import FWAction, FiretaskBase, Firework, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.fireworks.core import OptimizeFW
from mpinterfaces.transformations import get_aligned_lattices, \
    get_interface#generate_all_configs
from triboflow.utils.database import GetBulkFromDB, GetSlabFromDB, \
    GetHighLevelDB, GetDB
from triboflow.utils.vasp_tools import GetCustomVaspRelaxSettings
from triboflow.utils.structure_manipulation import InterfaceName


@explicit_serialize
class FT_StartSlabRelaxSWF(FiretaskBase):
    """Start a subworkflow as a detour to make a slab and relax it.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    slab_struct_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'unrelaxed_slab'.
    relax_type : str
        Relaxation type for the GetCustomVaspRelaxSettings helper_function.
    bulk_struct_name : str, optional
        Name of the bulk structure in the bulk database (material is
        identified by mp_id, but there might be differnt structures of the
        same material.) Defaults to 'structure_equiVol'.
    slab_out_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'relaxed_slab'.
        
    Returns
    -------
        Starts a new subworkflow as a detour to the current workflow.
    """
    
    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'slab_struct_name', 'relax_type',
                       'bulk_struct_name', 'slab_out_name']
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import MakeAndRelaxSlab_SWF
        mp_id = self.get('mp_id')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
            miller_str = self['miller']
        else:
            miller = self['miller']
            miller_str = ''.join(str(s) for s in self['miller'])
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        relax_type = self.get('relax_type', 'slab_pos_relax')
        bulk_name = self.get('bulk_struct_name', 'structure_equiVol')
        slab_out_name = self.get('slab_out_name', 'relaxed_slab')
        
        WF = MakeAndRelaxSlab_SWF(mp_id, miller, functional,
                                  relax_type = relax_type,
                                  bulk_struct_name = bulk_name,
                                  slab_struct_name = slab_name,
                                  out_struct_name = slab_out_name,
                                  spec = fw_spec)
        return FWAction(detours = WF)
    
@explicit_serialize
class FT_GetRelaxedSlab(FiretaskBase):
    """Get the relaxed structure from the DB, and put a Slab into the high-level DB.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    struct_out_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'relaxed_slab'.
    
    Returns
    -------
        Pymatgen Slab into the high-level DB.
    """
    
    required_params = ['mp_id', 'miller', 'functional', 'tag']
    optional_params = ['db_file', 'struct_out_name']
    def run_task(self, fw_spec):
        mp_id = self.get('mp_id')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
            miller_str = self['miller']
        else:
            miller = self['miller']
            miller_str = ''.join(str(s) for s in self['miller'])
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        out_name = self.get('struct_out_name', 'relaxed_slab')
        
        # Check if a relaxed slab is already in the DB entry
        slab_data = GetSlabFromDB(mp_id, db_file, miller, functional)
        if out_name not in slab_data:
            
            #Get results from OptimizeFW
            DB = GetDB(db_file)
            vasp_calc = DB.tasks.find_one({'task_label': self['tag']})
            relaxed_slab = Structure.from_dict(vasp_calc['output']['structure'])
        
            HL_DB = GetHighLevelDB(db_file)
            coll = HL_DB[functional+'.slab_data']
            
            slab_data = coll.find_one({'mpid': mp_id, 'miller': miller})
        
            #The following should be enough to transform the structure that
            #is coming out of the relaxation (relaxed_slab) into a Slab object.
                
            slab = Slab(relaxed_slab.lattice,
                        relaxed_slab.species_and_occu,
                        relaxed_slab.frac_coords,
                        miller,
                        Structure.from_sites(relaxed_slab, to_unit_cell=True),
                        shift=0,
                        scale_factor=[[1,0,0], [0,1,0], [0,0,1]],
                        site_properties=relaxed_slab.site_properties)
        
            coll.update_one({'mpid': mp_id, 'miller': miller},
                            {'$set': {out_name: slab.as_dict()}})
        
        return

            
@explicit_serialize
class FT_StartSlabRelax(FiretaskBase):
    """Fetch a Slab from the DB and relax it using an atomate OptimizeFW.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    slab_struct_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'unrelaxed_slab'.
    relax_type : str
        Relaxation type for the GetCustomVaspRelaxSettings helper_function.
    
    Returns
    -------
        Pymatgen Slab into the high-level DB.
    """
    
    required_params = ['mp_id', 'miller', 'functional', 'tag']
    optional_params = ['db_file', 'slab_struct_name', 'relax_type']
    def run_task(self, fw_spec):
        mp_id = self.get('mp_id')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
            miller_str = self['miller']
        else:
            miller = self['miller']
            miller_str = ''.join(str(s) for s in self['miller'])
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        tag = self.get('tag')
        relax_type = self.get('relax_type', 'slab_pos_relax')
        
        slab_data = GetSlabFromDB(mp_id, db_file, miller, functional)
        comp_params = slab_data.get('comp_parameters')
        slab_to_relax = Slab.from_dict(slab_data.get(slab_name))
        
        # Check if a relaxed slab is already in the DB entry
        if 'relaxed_slab' not in slab_data:
            vis = GetCustomVaspRelaxSettings(slab_to_relax, comp_params,
                                             relax_type)
        
            Optimize_FW = OptimizeFW(slab_to_relax, name=tag,
                                     vasp_input_set = vis,
                                     half_kpts_first_relax = True)
            wf_name = slab_data['formula']+miller_str+'_'+relax_type
            #Use add_modify_incar powerup to add KPAR and NCORE settings
            #based on env_chk in my_fworker.yaml
            Optimize_WF = add_modify_incar(Workflow([Optimize_FW],
                                                    name = wf_name))
            return FWAction(detours=Optimize_WF)

@explicit_serialize
class FT_MakeSlabInDB(FiretaskBase):
    """Makes a slab with certain orientation out of a bulk structure.

    The slab has a defined thickness and orientation.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    bulk_struct_name : str, optional
        Name of the bulk structure in the bulk database (material is
        identified by mp_id, but there might be differnt structures of the
        same material.) Defaults to 'structure_equiVol'.
    slab_struct_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'unrelaxed_slab'.
    
    Returns
    -------
        Pymatgen Slab into the high-level DB.
    """
    
    _fw_name = 'Make slab from bulk structure'
    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'bulk_struct_name', 'slab_struct_name']
    
    def run_task(self, fw_spec):
        
        mp_id = self.get('mp_id')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
        else:
            miller = self['miller']
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        bulk_name = self.get('bulk_struct_name', 'structure_equiVol')
        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        
        bulk_data = GetBulkFromDB(mp_id, db_file, functional)
        slab_data = GetSlabFromDB(mp_id, db_file, miller, functional)
        
        bulk_prim = Structure.from_dict(bulk_data[bulk_name])
        bulk_conv = SpacegroupAnalyzer(bulk_prim).get_conventional_standard_structure()
        
        min_thickness = slab_data.get('min_thickness', 10)
        min_vacuum = slab_data.get('min_min_vacuum', 25)
        
        SG = SlabGenerator(initial_structure = bulk_conv,
                           miller_index = miller,
                           center_slab = True,
                           primitive = True,
                           lll_reduce = True,
                           max_normal_search=max([abs(l) for l in miller]),
                           min_slab_size = min_thickness,
                           min_vacuum_size = min_vacuum)
        
# =============================================================================
# Here we need to do more Work! Now the first slab in the list is taken,
# which should be fine for monoatomic slabs, or if the methods only finds
# one slab which has no broken bonds, but in general, several slabs should be
# investigated for the lowest surface energy, and this one should be taken.
# =============================================================================
        slab = SG.get_slabs(bonds=None, ftol=0.1, tol=0.1, max_broken_bonds=0,
                            symmetrize=False, repair=False)[0]
        
        db = GetHighLevelDB(db_file)
        coll = db[functional+'.slab_data']
        
        slab_dict = monty.json.jsanitize(slab.as_dict(), allow_bson=True)
        
        coll.update_one({'mpid': mp_id, 'miller': miller},
                        {'$set': {'unrelaxed_slab': slab_dict}})
        
        
        return 



@explicit_serialize
class FT_MakeHeteroStructure(FiretaskBase):
    """Matches two slab systems to form a heterogeneous interface.

    If the match fails, the accuracy criteria are relaxed in steps of 5% until
    a match is found.
    
    Parameters
    ----------
    mp_id_1 : str
        ID number from the MP to identify material 1.
    miller_1 : list of int or str
        Miller indices for identification of the slab of material 1. Either
        single str e.g. '111', or a list of int e.g. [1, 1, 1].
    mp_id_2 : str
        ID number from the MP for material 2.
    miller_2 : list of int or str
        Miller indices for identification of the slab of material 2. Either
        single str e.g. '111', or a list of int e.g. [1, 1, 1].
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.

    Returns
    -------
    Matched interface structure and bottom and top slabs in interface_data DB.
    """
    
    _fw_name = 'Make Hetero Structure'
    required_params = ['mp_id_1', 'miller_1', 'mp_id_2', 'miller_2',
                       'functional']
    optional_params = ['db_file']
    
    def run_task(self, fw_spec):

        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        if type(self['miller_1']) == str:
            miller_1 = [int(k) for k in list(self['miller_1'])]
        else:
            miller_1 = self['miller_1']
        if type(self['miller_2']) == str:
            miller_2 = [int(k) for k in list(self['miller_2'])]
        else:
            miller_2 = self['miller_2']
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
            
        db = GetHighLevelDB(db_file)
        slab_col = db[functional+'.slab_data']
        inter_col = db[functional+'.interface_data']
        
        interface_name = InterfaceName(mp_id_1, miller_1, mp_id_2, miller_2)
        
        inter_data = inter_col.find_one({'name': interface_name})
        inter_params = inter_data['interface_parameters']
        comp_params = inter_data['comp_parameters']
        
        if not inter_data.get('unrelaxed_structure'):
           
            slab_1_dict = GetSlabFromDB(mp_id_1, db_file, miller_1, functional)
            slab_2_dict = GetSlabFromDB(mp_id_2, db_file, miller_2, functional)
            slab_1 = Slab.from_dict(slab_1_dict['relaxed_slab'])
            slab_2 = Slab.from_dict(slab_2_dict['relaxed_slab'])
            
# =============================================================================
# Running into crashes for max_angle_diff > ~1.5 for MPInterfaces 2020.6.19,
# at least for certain interfaces. A match is found, but than there is a 
# LinAlgError("Singular matrix") error in forming the matched slabs?
# The following lines ensures that max_angle_diff > 1.5. This is not a great
# solution obviously. I also changed the default in FT_CheckInterfaceParamDict
# =============================================================================
            if inter_params['max_angle_diff'] > 1.5:
                inter_params['max_angle_diff'] = 1.5
        
            bottom_aligned, top_aligned = get_aligned_lattices(
                slab_1,
                slab_2,
                max_area = inter_params['max_area'],
                max_mismatch = inter_params['max_mismatch'],
                max_angle_diff = inter_params['max_angle_diff'],
                r1r2_tol = inter_params['r1r2_tol'])
        
            if bottom_aligned:
# ============================================================================
# TODO: Find out if this is actually useful for the PES to return
#       all the interface structures we need, or if another stacking function
#       should be written with another input of lateral shifts.
#       Check out AdsorbateSiteFinder of pymatgen.analysis.adsorption!!
#       If generate all configs will work, we just have to remove the
#       [0] index and move the whole list of structures to the spec.
# EDIT (19.08.2020): It seems as if the code of MPInterfaces has been changed.
#       Now the relevant function is called get_interface and not
#       generate_all_configs any more, and only one Interface is returned... 
# ============================================================================
# ============================================================================
                #hetero_interfaces = generate_all_configs(top_aligned,
                hetero_interfaces = get_interface(bottom_aligned,
                                                  top_aligned,
                                                  nlayers_2d = 1,
                                                  nlayers_substrate = 1,
                            separation = inter_params['interface_distance'])
        
                #inter_slab = hetero_interfaces[0]
                inter_slab = hetero_interfaces
                
                inter_dict = inter_slab.as_dict()
                bottom_dict = bottom_aligned.as_dict()
                top_dict = top_aligned.as_dict()
                inter_col.update_one({'name': interface_name},
                                 {'$set': {'unrelaxed_structure': inter_dict,
                                           'bottom_aligned': bottom_dict,
                                           'top_aligned': top_dict}})
            
            else:
                f = 1.05 # factor with wich to multiply the missmatch criteria
            
                new_params = {'max_area': inter_params['max_area'] * f,
                    'max_mismatch': inter_params['max_mismatch'] * f,
                    'max_angle_diff': inter_params['max_angle_diff'] * f,
                    'r1r2_tol': inter_params['r1r2_tol'] * f,
                    'interface_distance': inter_params['interface_distance']}
            
                if not inter_data.get('original_interface_params'):
                    inter_col.update_one({'name': interface_name},
                        {'$set': {'original_interface_params': inter_params}})
                    
                inter_col.update_one({'name': interface_name},
                        {'$set': {'interface_parameters': new_params}})
            
            
            
                new_fw = Firework(FT_MakeHeteroStructure(mp_id_1 = mp_id_1,
                                                         miller_1 = miller_1,
                                                         mp_id_2 = mp_id_2,
                                                         miller_2 = miller_2,
                                                         db_file = db_file))
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



