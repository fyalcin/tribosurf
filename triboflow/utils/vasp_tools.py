import os
import subprocess

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet, MPStaticSet

from triboflow.utils.file_manipulation import RemoveMatchingFiles

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
    uis['NELMIN'] = 5
    uis['EDIFF'] = 0.5E-5
    uis['LAECHG'] = '.FALSE.'
    
    if structure.num_sites < 20:
        uis['LREAL'] = '.FALSE.'
    
    if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
        uis['NELMDL'] = -15
        uis['EDIFFG'] = -0.02
        uis['NELM'] = 200
        # Turn on linear mixing
        # uis['AMIX'] = 0.2
        # uis['BMIX'] = 0.0001
        # uis['AMIX_MAG'] = 0.8
        # uis['BMIX_MAG'] = 0.0001
        
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