"""A collection of HelperFunctions to be used for the FireFlow project."""

import os, glob
import subprocess
import pymongo
import io
import numpy as np
from PIL import Image
from pymatgen import MPRester
from pymatgen.core.surface import Slab
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.inputs import Kpoints, Poscar
from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet, MPStaticSet
from atomate.vasp.database import VaspCalcDb


def SlabFromStructure(miller, structure):
    """Returns a pymatgen.core.surface.Slab from a pymatgen structure.

    Parameters
    ----------
    miller : list of int
        Miller indices given as a list of integers
    structure : pymatgen.core.structure.Structure
        Structure to be converted into a Slab

    Returns
    -------
    pymatgen.core.surface.Slab
        The input structure converted to a Slab

    """
    return Slab(lattice = structure.lattice,
                species = structure.species_and_occu,
                coords = structure.frac_coords,
                miller_index = miller,
                oriented_unit_cell = structure,
                shift = 0,
                scale_factor = np.eye(3, dtype=np.int),
                site_properties = structure.site_properties)

def ConvertBytesToImage(bytes_object):
    """Convert an image saved as bytes for storing in MongoDB to an PIL image.

    Parameters
    ----------
    bytes_object : bytes
        Image converted to bytes for storage in a MongoDB database.

    Returns
    -------
    pil_img : PIL.PngImagePlugin.PngImageFile
        Image in python image library format ready to be viewed or saved.

    """
    pil_img = Image.open(io.BytesIO(bytes_object))
    return pil_img

def ConvertFigToBytes(path_to_fig):
    """Convert an image to bytes for starage in MongoDB database.

    Parameters
    ----------
    path_to_fig : Str
        Path to the input figure.

    Returns
    -------
    bytes
        image encoded in bytes ready for storage in MongoDB.

    """
    im = Image.open(path_to_fig)
    image_bytes = io.BytesIO()
    im.save(image_bytes, format='png')
    return image_bytes.getvalue()

def CleanUpSitePorperties(structure):
    """
    Cleans up site_properties of structures that contain NoneTypes.
    
    If an interface is created from two different structures, it is possible
    that some site properties like magmom are not set for both structures.
    This can lead later to problems since they are replaced by None.
    This function replaces NoneTypes with 0.0 for magmom and deletes all other
    site_properties if None entries are found in it.


    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Input structure

    Returns
    -------
    struct : pymatgen.core.structure.Structure
        Output structure

    """
    struct = structure.copy()
    for key in struct.site_properties.keys():
        if key == 'magmom':
            new_magmom = []
            for m in struct.site_properties[key]:
                if m == None:
                    new_magmom.append(0.0)
                else:
                    new_magmom.append(m)
            struct.add_site_property('magmom', new_magmom)
        else:
            struct.remove_site_property(key)
    return struct

def StackAlignedSlabs(bottom_slab, top_slab, top_shift=[0,0,0]):
    """
    Combine slabs that are centered around 0 into a single structure.
    
    Optionally shift the top slab by a vector of cartesian coordinates.

    Parameters
    ----------
    bottom_slab : pymatgen.core.structure.Structure
        Bottom slab.
    top_slab : pymatgen.core.structure.Structure
        Top slab.
    top_shift : list of 3 floats, optional
        Vector of caresian coordinates with which to shift the top slab.
        The default is [0,0,0].

    Returns
    -------
    interface : pymatgen.core.structure.Structure
        DESCRIPTION.

    """
    interface = bottom_slab.copy()
    t_copy = top_slab.copy()
    
    t_copy.translate_sites(indices=range(len(t_copy.sites)),
                           vector=top_shift,
                           frac_coords=False, to_unit_cell=False)
    
    for s in t_copy.sites:
        new_site = PeriodicSite(lattice=interface.lattice,
                                coords=s.frac_coords,
                                coords_are_cartesian=False,
                                species=s.species)
        interface.sites.append(new_site)
    
    return interface

def ReCenterAlignedSlabs(top_slab, bottom_slab, d=2.5):
    """
    Center two slabs around z=0 and give them the distance d.

    Parameters
    ----------
    top_slab : pymatgen.core.structure.Structure
        The slab that should be on top.
    bottom_slab : pymatgen.core.structure.Structure
        The slab that should be on the bottom.
    d : float, optional
        The desired distance between the slabs. The default is 2.5.

    Returns
    -------
    t_copy : pymatgen.core.structure.Structure
        Top slab that is shifted so that the lowest atom is at +d/2
    b_copy : pymatgen.core.structure.Structure
        Bottom slab that is shifted so that the topmost atom is at -d/2

    """
    t_copy = top_slab.copy()
    b_copy = bottom_slab.copy()
    top_zs=[]
    bot_zs=[]
    for s in t_copy.sites:
        top_zs.append(s.coords[-1])
    top_shift = -min(top_zs) + d/2
    
    for s in b_copy.sites:
        bot_zs.append(s.coords[-1])
    bot_shift = -max(bot_zs) - d/2

    t_copy.translate_sites(indices=range(len(t_copy.sites)),
                           vector=[0,0,top_shift],
                           frac_coords=False, to_unit_cell=False)
    b_copy.translate_sites(indices=range(len(b_copy.sites)),
                           vector=[0,0,bot_shift],
                           frac_coords=False, to_unit_cell=False)
    return t_copy, b_copy

def GetEmin(potcar):
    """
    Return the minimal recommended  energy cutoff for a given Potcar object.
    
    Unfortunately I don not think that pymatgen Potcars can access the data
    for the recommended cutoffs directly, so I write out a file first and
    then scan for the ENMAX lines. The largest EMIN value found is returned
    in the end.

    Parameters
    ----------
    potcar : pymatgen.io.vasp.inputs.Potcar
        The pymatgen representation of a VASP POTCAR file.

    Returns
    -------
    float
        The largest EMIN value of all species present in the potcar.

    """
    potcar.write_file('temp_potcar')
    with open('temp_potcar', 'r') as pot:
        emin = []
        for l in pot:
            if l.strip().startswith('ENMAX'):
                emin.append(float(l.split()[-2]))
    os.remove('temp_potcar')
    return max(emin)
                    

def RemoveMatchingFiles(list_of_patterns):
    """
    Remove all files matching the patterns (wildcards) in the current
    directory.
    """
    remove_list=[]
    for pattern in list_of_patterns:
        remove_list.extend(glob.glob(pattern))
    if remove_list is []:
        return
    else:
        for file_to_remove in remove_list:
            os.remove(file_to_remove)
        return

def GetGeneralizedKmesh(structure, k_dist, RemoveSymm=False):
    """Get a generalized Monkhorst pack mesh for a given structure.
    
    Prepares the necessary files (POSCAR, PRECALC, and, if the structure has
    a 'magmom' site property also INCAR) for the K-Point Grid Generator of the
    Mueller group at John Hopkins http://muellergroup.jhu.edu/K-Points.html
    Runs the getKPoints script and reads the KPOINT file produced into a
    pymatgen Kpoints object.
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Structure for which the kpoints grid is to be generated.
    k_dist : float
        The minimum allowed distance between lattice points on the real-space
        superlattice. This determines the density of the k-point grid. A larger
        value will result in denser grids.

    Returns
    -------
    KPTS : pymatgen.io.vasp.inputs.Kpoints
        Pymatgen Kpoints object representing a generalized MP mesh.

    """
    precalc = ['MINDISTANCE = '+str(k_dist),
               'MINTOTALKPOINTS = 4',
               'GAPDISTANCE = 6 ',
               'MONOCLINIC_SEARCH_DEPTH = 2500',
               'TRICLINIC_SEARCH_DEPTH = 1500']
    if RemoveSymm in ['STRUCTURAL', 'TIME_REVERSAL', 'ALL']:
        precalc.append('REMOVE_SYMMETRY = '+RemoveSymm)
    with open('PRECALC', 'w') as out:
        for line in precalc:
            out.write(line+'\n')
            
    magmom_list = structure.site_properties.get('magmom')
    if magmom_list:
        with open('INCAR', 'w') as out:
            out.write('ISPIN = 2')
            out.write('MAGMOM = '+' '.join(str(m) for m in magmom_list))
    
    structure.to(fmt='poscar', filename='POSCAR')
    get_kpoints_file = subprocess.Popen('getKPoints')
    get_kpoints_file.communicate()
    KPTS = Kpoints().from_file('KPOINTS')
    RemoveMatchingFiles(['KPOINTS*', 'POSCAR*', 'INCAR', 'PRECALC'])
    return KPTS

def InterfaceName(mp_id_1, miller_1, mp_id_2, miller_2):
    """Return a name for an interface based on MP-IDs and miller indices.

    Parameters
    ----------
    mp_id_1 : str
        MP-ID of the first material.
    miller_1 : list of int, or str
        Miller indices of the first material given either as list or str with
        3 letters.
    mp_id_2 : str
        MP-ID of the second material.
    miller_2 : list of int, or str
        Miller indices of the second material given either as list or str with
        3 letters.

    Returns
    -------
    name : str
        Unique name for the interface of two slabs.

    """
    f1 = GetPropertyFromMP(mp_id_1, 'pretty_formula')
    f2 = GetPropertyFromMP(mp_id_2, 'pretty_formula')
    if type(miller_1) is list:
        m1 = ''.join(str(s) for s in miller_1)
    else:
        m1 = miller_1
    if type(miller_2) is list:
        m2 = ''.join(str(s) for s in miller_2)
    else:
        m2 = miller_2
    n1 = min(f1+m1, f2+m2)
    n2 = max(f1+m1, f2+m2)
    ids = min(mp_id_1+'_'+mp_id_2, mp_id_2+'_'+mp_id_1)
    name = '_'.join((n1, n2, ids))
    return name

def GetHighLevelDB(db_file):
    """Return the triboflow MongoDB database.
    
    Parameters
    ----------
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database.
        The default is '/home/mwo/FireWorks/config/db.json'.
    
    Returns
    -------
    db : pymongo.database.Database
        Triboflow high-level MongoDB database accessed by pymongo.
    """
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    
    return tribo_db
    
def AddBulkToDB(structure, mp_id, db_file, functional):
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.bulk_data'
    col = tribo_db[col_name]
    formula = structure.composition.reduced_formula
    
    if col.find_one({'mpid': mp_id}):
        print('{} bulk can not be added to bulk_data collection because a '
              'material with MP-ID {} is already present in the {} '
              'collection!'.format(formula, mp_id, functional))
        return
    
    col.insert_one({'mpid': mp_id,
                    'formula': formula,
                    'structure_fromMP': structure.as_dict()})
    return

def GetBulkFromDB(mp_id, db_file, functional):
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.bulk_data'
    col = tribo_db[col_name]
    out_dict = col.find_one({'mpid': mp_id})
    if out_dict:
        return out_dict
    else:
        raise IOError('No bulk material with MP-ID {} is found in the'
                      '{}.bulk_data collection.'.format(mp_id, functional))
        
def GetSlabFromDB(mp_id, db_file, miller, functional):
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.slab_data'
    col = tribo_db[col_name]
    out_dict = col.find_one({'mpid': mp_id, 'miller': miller})
    if out_dict:
        return out_dict
    else:
        raise IOError('No slab with MP-ID {} and miller indices {} was found'
                      ' in the {}.slab_data collection.'
                      .format(mp_id, miller, functional))
        
def GetInterfaceFromDB(name, db_file, functional):
    db = GetDB(db_file)
    tribo_db = db.client.triboflow
    col_name=functional+'.interface_data'
    col = tribo_db[col_name]
    out_dict = col.find_one({'name': name})
    if out_dict:
        return out_dict
    else:
        raise IOError('No interface with name {} was found in the '
                      '{}.interface_data collection.'.format(name, functional))
        

def GetLastBMDatafromDB(formula, db_file='/home/mwo/FireWorks/config/db.json'):
    """Query the FireWorks database for the last bulk modulus data for formula.
    
    A query is made to the eos collection in the Fireworks database for a
    given chemical formula. The results are sorted ascending by creation date
    and than the last one is returned as a dictionary.

    Parameters
    ----------
    formula : str
        Chemical formula for which the query is made.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database.
        The default is '/home/mwo/FireWorks/config/db.json'.

    Returns
    -------
    dict
        Data generated by the get_wf_bulk_modulus workflow of atomate.

    """
    db = GetDB(db_file)
    cursor = db.eos.find({'formula_pretty': formula}).sort(
        'created_at', pymongo.DESCENDING)
    return cursor[0]

def GetDB(db_file='/home/mwo/FireWorks/config/db.json'):
    """Connect to the MongoDB database specified in the db_file.
    
    Parameters
    ----------
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database.
        The default is '/home/mwo/FireWorks/config/db.json'.

    Returns
    -------
    db : pymongo.database.Database
        MongoDB database accessed by pymongo

    """
    VaspCalcDB = VaspCalcDb.from_db_file(db_file)
    return VaspCalcDB.db

def IsListConverged(input_list, tol, n=3):
    """Check if the last n values of an array are within tol of each other.
    

    Parameters
    ----------
    input_list : list of float
        Total energies to be checked for convergence
    tol : float
        Tolerance for the convergence.
    n : int, optional
        Number of entries at the end of energy_list that have to be within
        etol for the list to be considered converged. The default is 3.

    Returns
    -------
    Bool
        True if input_list is converged, False otherwise.

    """
    if len(input_list) <= n:
        return False
    else:
        check_list = [False]*n
        l = input_list.copy()
        l.reverse()
        for i, b in enumerate(check_list):
            if abs(l[0]-l[i+1]) < tol:
                check_list[i] = True
        return all(check_list)

def GetCustomVaspStaticSettings(structure, comp_parameters, static_type):
    """Make custom vasp settings for static calculations.
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Structure to be treated
    comp_parameters : dict
        Computational parameters dictionary which is usually created partly
        by the user and then filled automatically with defaults and after
        convergence tests.
    statc_type : str
        Specifies what system is treated in what way. Check 'allowed_types'
        for a list of choices.

    Raises
    ------
    SystemExit
        If a non-supported static_type is passed, the process terminates.

    Returns
    -------
    vis : pymatgen.io.vasp.sets.MPStaticSet
        
    """

    allowed_types = ['bulk_from_scratch', 'bulk_follow_up', 'bulk_nscf',
                     'slab_from_scratch', 'slab_follow_up', 'slab_nscf']
    
    if static_type not in allowed_types:
        raise SystemExit('static type is not known. Please select from: {}'
                         .format(allowed_types))
        
    #Set user incar settings:
    uis = {}
    uis['NEDOS'] = 3001
    uis['PREC'] = 'Accurate'
    uis['GGA_COMPAT'] = '.FALSE.'
    uis['LASPH'] = '.TRUE.'
    uis['LORBIT'] = 11
    uis['NELMIN'] = 4
    uis['SIGMA'] = 0.05
    uis['ISMEAR'] = -5
    uis['EDIFF'] = 1.0e-6
    
    if static_type.endswith('from_scratch'):
        uis['ICHARG'] = 2
        uis['LAECHG'] = '.FALSE.'
    
    if structure.num_sites < 20:
        uis['LREAL'] = '.FALSE.'
        
    if comp_parameters.get('functional') == 'SCAN':
        uis['ISMEAR'] = 0
        uis['SIGMA'] = 0.1
    
    if static_type.startswith('slab_'):
        uis['NELMDL'] = -15
        uis['NELM'] = 200
    elif comp_parameters.get('functional') == 'SCAN':
        uis['NELMDL'] = -10
    else:
        uis['NELMDL'] = -6
        
    if 'encut' in comp_parameters:
        uis['ENCUT'] = comp_parameters['encut']
        
    if 'use_spin' in comp_parameters:
        if comp_parameters['use_spin']:
            uis['ISPIN'] = 2
        else:
            uis['ISPIN'] = 1
        
    #set van der Waals functional. Note that as of now, 'functional' must be
    #specified for vdw to work!
    if set(('use_vdw', 'functional')) <= comp_parameters.keys():
        if comp_parameters['use_vdw']:
            if comp_parameters['functional'] == 'SCAN':
                vdw = 'rVV10'
                uis['LUSE_VDW'] = '.TRUE.'
                uis['BPARAM'] = 15.7
            else:
                vdw = 'optB86b'
        else:
            vdw = None
    else:
        vdw = None
    
    if comp_parameters.get('functional') == 'SCAN':
        uis['METAGGA'] = 'SCAN'
        uis['ALGO'] = 'All'
        
    if static_type.endswith('follow_up'):
        uis['ISTART'] = 1
        uis['LREAL'] = '.FALSE.'
        uis['NELMDL'] = -1
    elif static_type.endswith('nsfc'):
        uis['ISTART'] = 1
        uis['LREAL'] = '.FALSE.'
        uis['ICHARG'] = 11
        uis['NELMDL'] = -1
        
    if 'k_dens' in comp_parameters:
        kpoints = Kpoints.automatic_gamma_density(structure,
                                              comp_parameters['k_dens'])
    else:
        #if no k-density is supplied in the comp_parameters, use a large value
        kpoints = Kpoints.automatic_gamma_density(structure, 5000)
    
    #If a structure has a vacuum layer, set the third kpoints devision to 1
    #by force.
    if static_type.startswith('slab_'):
        uks = Kpoints(comment=kpoints.comment+'  adjusted for slabs',
                      num_kpts=0,
                      kpts=[[kpoints.kpts[0][0], kpoints.kpts[0][1], 1]])
    else:
        uks = kpoints
        
    # Set vasp input set (currently none available for static SCAN!)
    vis = MPStaticSet(structure, user_incar_settings = uis, vdw = vdw,
                      user_kpoints_settings = uks,
                      user_potcar_functional = 'PBE_54')
        
    return vis

def GetCustomVaspRelaxSettings(structure, comp_parameters, relax_type):
    """Make custom vasp settings for relaxations.
    
    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Structure to be relaxed.
    comp_parameters : dict
        Computational parameters dictionary which is usually created partly
        by the user and then filled automatically with defaults and after
        convergence tests.
    relax_type : str
        Specifies what is to be relaxed in what way. Check 'allowed_types'
        for a list of choices.

    Raises
    ------
    SystemExit
        If a non-supported relax_type is passed, the process terminates.

    Returns
    -------
    vis : pymatgen.io.vasp.sets.MPScanRelaxSet or
          pymatgen.io.vasp.sets.MPRelaxSet, depending on input
        A vasp input set for pymatgen.

    """
    allowed_types = ['bulk_full_relax', 'bulk_vol_relax', 'bulk_pos_relax',
                     'bulk_shape_relax',
                     'slab_shape_relax', 'slab_pos_relax',
                     'interface_shape_relax', 'interface_pos_relax',
                     'interface_z_relax']
    
    if relax_type not in allowed_types:
        raise SystemExit('relax type is not known. Please select from: {}'
                         .format(allowed_types))
    
    #Set user incar settings:
    uis = {}
    uis['NEDOS'] = 3001
    uis['PREC'] = 'Accurate'
    uis['GGA_COMPAT'] = '.FALSE.'
    uis['LASPH'] = '.TRUE.'
    uis['LORBIT'] = 11
    uis['MAXMIX'] = 100
    uis['NELMIN'] = 4
    uis['EDIFF'] = 1.0E-5
    uis['LAECHG'] = '.FALSE.'
    
    if structure.num_sites < 20:
        uis['LREAL'] = '.FALSE.'
    
    if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
        uis['NELMDL'] = -15
        uis['EDIFFG'] = -0.02
        uis['NELM'] = 200
    else:
        uis['NELMDL'] = -6
        uis['EDIFFG'] = -0.01
    
    if relax_type.endswith('full_relax'):
        uis['ISIF'] = 3
    elif relax_type.endswith('pos_relax'):
        uis['ISIF'] = 2
    elif relax_type.endswith('z_relax'):
        uis['ISIF'] = 2
        #Set up selective dynamics array for the structrues site property
        sd_array = []
        for i in range(len(structure.sites)):
            sd_array.append([False, False, True])
        structure.add_site_property('selective_dynamics', sd_array)
    elif relax_type.endswith('vol_relax'):
        uis['ISIF'] = 7
    elif relax_type.endswith('shape_relax'):
        uis['ISIF'] = 5
    
    if 'encut' in comp_parameters:
        uis['ENCUT'] = comp_parameters['encut']
        
    if 'use_spin' in comp_parameters:
        if comp_parameters['use_spin']:
            uis['ISPIN'] = 2
        else:
            uis['ISPIN'] = 1
    
    if 'is_metal' in comp_parameters:
        if comp_parameters['is_metal']:
            uis['SIGMA'] = 0.2
            uis['ISMEAR'] = 2
        else:
            uis['SIGMA'] = 0.05
            uis['ISMEAR'] = -5
    else:
        uis['SIGMA'] = 0.1
        uis['ISMEAR'] = 0
        
        
    #set van der Waals functional. Note that as of now, 'functional' must be
    #specified for vdw to work!
    if set(('use_vdw', 'functional')) <= comp_parameters.keys():
        if comp_parameters['use_vdw']:
            if comp_parameters['functional'] == 'SCAN':
                vdw = 'rVV10'
            else:
                vdw = 'optB86b'
        else:
            vdw = None
    else:
        vdw = None
        
    if 'k_dens' in comp_parameters:
        kpoints = Kpoints.automatic_gamma_density(structure,
                                              comp_parameters['k_dens'])
    else:
        #if no k-density is supplied in the comp_parameters, use a large value
        kpoints = Kpoints.automatic_gamma_density(structure, 5000)

    #If a structure has a vacuum layer, set the third kpoints devision to 1
    #by force.
    if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
        uks = Kpoints(comment=kpoints.comment+'  adjusted for slabs',
                      num_kpts=0,
                      kpts=[[kpoints.kpts[0][0], kpoints.kpts[0][1], 1]])
    else:
        uks = kpoints
        
    if 'functional' in comp_parameters:
        if comp_parameters['functional'] == 'SCAN':
            #Algo All does not play well with tetrahedron method
            if 'is_metal' in comp_parameters:
                if not comp_parameters['is_metal']:
                    uis['SIGMA'] = 0.1
                    uis['ISMEAR'] = 0
            vis = MPScanRelaxSet(structure, user_incar_settings = uis,
                                 vdw = vdw, user_kpoints_settings = uks,
                                 user_potcar_functional = 'PBE_54')
        else:
            vis = MPRelaxSet(structure, user_incar_settings = uis, vdw = vdw,
                             user_kpoints_settings = uks,
                             user_potcar_functional = 'PBE_54')
    else:
       vis = MPRelaxSet(structure, user_incar_settings = uis, vdw = vdw,
                        user_kpoints_settings = uks,
                        user_potcar_functional = 'PBE_54')
        
    return vis


def GetPropertyFromMP(MP_ID, prop):
    """Get a certain property for a single material from the MP database.
    
    Return the selected property of the selected material from the Material
    Projects database using the MPRester from Pymatgen.

    Parameters
    ----------
    MP_ID : str
        Valid materials project ID.
    prop : str
        Property for which to query. If not in the supported list, a warning
        will be printed.

    Returns
    -------
    variable type
        The query result from the MP database. Type depends on the query.
        If nothing is found or the property is not in the list, NoneType is
        returned.
    """
    with MPRester() as mpr:
        
        if prop not in mpr.supported_task_properties:
            print('{} is not in the list of supported task properties!\n'
                  'This is the supported list:\n{}'
                  .format(prop, mpr.supported_task_properties))
            return None
        else:
            query=mpr.query(criteria={'material_id': MP_ID},
                            properties=[prop])
            return query[0][prop]

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
