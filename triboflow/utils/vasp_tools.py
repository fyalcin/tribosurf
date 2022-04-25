import os
import subprocess
import math

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import (MPRelaxSet, MPScanRelaxSet, MPStaticSet,
                                   MPScanStaticSet)

from triboflow.utils.file_manipulation import remove_matching_files


class MeshFromDensity:
    """
    Class to find classic Monkhorst-Pack meshes which may
    or my not be Gamma-Centred from a given k-point density.
    Provides also capabilites to check if meshes generated from similar
    densities are equivalent. Meshes for slabs are adapted according to
    the lower symmetry conditions.
    
    """

    def __init__(self,
                 structure,
                 target_density,
                 compare_density=10.0,
                 is_slab='Auto',
                 min_vac=6.0,
                 force_gamma=False):
        """Initializes the class and sets internal variables.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            The pymatgen representation of a structure
        target_density : float
            Desired minimal density of kpoints along each reciprocal lattice
            vector in 1/Angstrom.
        compare_density : float
            Density for which the mesh of the target density is to be compared
            to. in 1/Angstom. The default is 10.0
        is_slab : bool/str
            If the passed structure is to be considered to be a slab. If "True"
            it is considered to be a slab, if "False" it is considered to be
            bulk, if "Auto", it is attempted to find out. The default is
            "Auto".
        min_vac : float
            Minimum vacuum distance for the slab detection if is_slab is set to
            "Auto". The default is 6.0
        force_gamma : bool
            If Gamma-centred meshes are to be forces (not applicable for
            generalized meshes!). The default is False.

        """

        self.struct = structure.copy()
        self.dens = target_density
        self.compare_dens = compare_density
        self.min_vac = min_vac
        self.force_gamma = force_gamma
        self.klm = structure.lattice.reciprocal_lattice.abc

        if is_slab == True:
            self.slab = True
        elif is_slab in ['Auto', 'auto', 'automatic', 'Automatic']:
            self.slab = 'detect_automatically'
        else:
            self.slab = False

    def __make_mesh(self, density):
        """Return the subdivisions along each lattice vector.
        
        Consider also if the structure is a slab.
        

        Parameters
        ----------
        density : float
            Desired minimal density of kpoints along each reciprocal lattice
            vector in 1/Angstrom.

        Returns
        -------
        k1 : int
            Kpoint devisions along b1
        k2 : int
            Kpoint devisions along b2
        k3 : int
            Kpoint devisions along b3

        """

        k, l, m = self.klm
        k1 = math.ceil(k * density)
        k2 = math.ceil(l * density)
        k3 = math.ceil(m * density)
        if self._is_slab():
            k3 = 1
        return (k1, k2, k3)

    def _is_slab(self):
        """Figures out if the passed structure is a slab.
        
        Automatic detection might fail for slabs that are set up in a none
        standard way!
        
        Returns
        -------
        bool
            True if the structure is considered a slab, False if not.

        """

        if self.slab == True:
            return True
        elif self.slab == 'detect_automatically':
            z_axis = self.struct.lattice.c
            z_coords = []
            for s in self.struct.sites:
                z_coords.append(s.coords[2])

            thickness = max(z_coords) - min(z_coords)
            if z_axis - thickness >= self.min_vac:
                return True
            else:
                return False
        else:
            return False

    def get_kpoints(self):
        """Return a Kpoint object with the desired density of kpoints.

        Returns
        -------
        kpoints : pymatgen.io.vasp.Kpoints
            Monkhorst-Pack or Gamma centered mesh.

        """

        mesh = self.__make_mesh(self.dens)
        is_hexagonal = self.struct.lattice.is_hexagonal()
        # has_odd = any(i % 2 == 1 for i in mesh)

        if is_hexagonal or self.force_gamma:
            kpoints = Kpoints.gamma_automatic(kpts=mesh)
        else:
            kpoints = Kpoints.monkhorst_automatic(kpts=mesh)

        return kpoints

    def are_meshes_the_same(self):
        """Compares conventional Monkhorst-Pack meshes and Gamma centered meshes.
        
        To test if a different target density actually provides a different
        mesh than a reference density.

        Returns
        -------
        bool
            True if meshes are the same, False otherwise.

        """

        mesh_1 = self.__make_mesh(self.dens)
        mesh_2 = self.__make_mesh(self.compare_dens)
        if mesh_1 == mesh_2:
            return True
        else:
            return False


def get_emin_and_emax(potcar):
    """
    Return ENMIN and ENMAX energy cutoff for a given Potcar object.
    
    Unfortunately I don not think that pymatgen Potcars can access the data
    for the recommended cutoffs directly, so I write out a file first and
    then scan for the ENMAX lines. The largest EMIN and ENMAX values found are
    returned in the end.

    Parameters
    ----------
    potcar : pymatgen.io.vasp.inputs.Potcar
        The pymatgen representation of a VASP POTCAR file.

    Returns
    -------
    dict
        The largest EMIN and ENMAX value of all species present in the potcar.

    """

    potcar.write_file('temp_potcar')
    with open('temp_potcar', 'r') as pot:
        emin = []
        emax = []
        for l in pot:
            if l.strip().startswith('ENMAX'):
                emin.append(float(l.split()[-2]))
                emax.append(float(l.split()[2][:-1]))
    os.remove('temp_potcar')
    return {'ENMIN': max(emin), 'ENMAX': max(emax)}


def get_generalized_kmesh(structure, k_dist, RemoveSymm=False, Vasp6=True):
    """Get a generalized Monkhorst Pack mesh for a given structure.
    
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

    precalc = ['MINDISTANCE = ' + str(k_dist),
               'MINTOTALKPOINTS = 4',
               'GAPDISTANCE = 6 ',
               'MONOCLINIC_SEARCH_DEPTH = 2500',
               'TRICLINIC_SEARCH_DEPTH = 1500']
    if Vasp6:
        precalc.append('WRITE_LATTICE_VECTORS = True')

    if RemoveSymm in ['STRUCTURAL', 'TIME_REVERSAL', 'ALL']:
        precalc.append('REMOVE_SYMMETRY = ' + RemoveSymm)
    with open('PRECALC', 'w') as out:
        for line in precalc:
            out.write(line + '\n')

    magmom_list = structure.site_properties.get('magmom')
    if magmom_list:
        with open('INCAR', 'w') as out:
            out.write('ISPIN = 2')
            out.write('MAGMOM = ' + ' '.join(str(m) for m in magmom_list))

    structure.to(fmt='poscar', filename='POSCAR')
    get_kpoints_file = subprocess.Popen('getKPoints')
    get_kpoints_file.communicate()
    KPTS = Kpoints().from_file('KPOINTS')
    remove_matching_files(['KPOINTS*', 'POSCAR*', 'INCAR', 'PRECALC'])

    return KPTS


def get_custom_vasp_static_settings(structure, comp_parameters, static_type,
                                    k_dens_default=8.5, ups={'W': 'W_sv'}):
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
    k_dens_default : float, optional
        Specifies the default kpoint density if no k_dens or kspacing key
        is found in the comp_parameters dictionary. The default is 12.5
    ups : dict, optional
        Specify user potcar settings, e.g. which potcar to choose for which
        element. The default is there to fix an issue with tungsten, where
        MP uses W_pv, which is depreciated by VASP and replaced with W_sv.
        The default is {'W': 'W_sv'}
        
    Raises
    ------
    SystemExit
        If a non-supported static_type is passed, the process terminates.

    Returns
    -------
    vis : pymatgen.io.vasp.sets.MPStaticSet
        
    """

    allowed_types = ['bulk_from_scratch', 'bulk_follow_up', 'bulk_nscf',
                     'slab_from_scratch', 'slab_follow_up', 'slab_nscf',
                     'bulk_epsilon_from_scratch', 'bulk_epsilon_follow_up']

    if static_type not in allowed_types:
        raise SystemExit('static type is not known. Please select from: {}'
                         .format(allowed_types))

    SCAN_list = ['scan', 'rscan', 'r2scan', 'Scan', 'Rrscan', 'R2scan',
                 'SCAN', 'RSCAN', 'R2SCAN']

    # Set user incar settings:
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
    uis['SYMPREC'] = 1e-04  # compat some issues that VASP 6.2 has with kpoint lattices
    #include electrons with l-quantum number up to 6 into the mixer. Helps with convergence
    #maybe include functionality that sets LMAXMIX dependent on periodic table group
    # e.g.:
    # for element in structure.composition.elements:
    #     if element.group == 3:
    #         uis['LMAXMIX'] = 6
    #     elif element.group in [4,5,6,7,8,9,10,11,12]:
    #         uis['LMAXMIX'] = 4
    uis['LMAXMIX'] = 6 

    if static_type.startswith('bulk_'):
        uis['ALGO'] = 'Fast'

    if static_type.endswith('from_scratch'):
        uis['ICHARG'] = 2
        uis['LAECHG'] = '.FALSE.'

    if structure.num_sites < 20:
        uis['LREAL'] = '.FALSE.'

    # Adjust mixing for slabs that have a very large c axis:
    if structure.lattice.matrix[-1, 1] > 50.0:
        uis['AMIN'] = 0.05

    if comp_parameters.get('functional') in SCAN_list:
        uis['ISMEAR'] = 0
        uis['SIGMA'] = 0.1

    if static_type.startswith('slab_'):
        uis['NELMDL'] = -15
        uis['NELM'] = 200
        # set dipole corrections.
        try:
            uis['DIPOL'] = list(structure.center_of_mass)
            uis['IDIPOL'] = 3
            uis['EPSILON'] = comp_parameters.get('epsilon') or 1.0
        except:
            uis['IDIPOL'] = 3
            uis['EPSILON'] = comp_parameters.get('epsilon') or 1.0
    elif comp_parameters.get('functional') in SCAN_list:
        uis['NELMDL'] = -10
    else:
        uis['NELMDL'] = -6
        uis['NELM'] = 200

    if 'encut' in comp_parameters:
        uis['ENCUT'] = comp_parameters['encut']

    if 'use_spin' in comp_parameters:
        if comp_parameters['use_spin']:
            uis['ISPIN'] = 2
        else:
            uis['ISPIN'] = 1

    # set van der Waals functional. Note that as of now, 'functional' must be
    # specified for vdw to work!
    if set(('use_vdw', 'functional')) <= comp_parameters.keys():
        if comp_parameters['use_vdw']:
            if comp_parameters.get('functional') in SCAN_list:
                vdw = 'rVV10'
            else:
                vdw = 'optB86b'
        else:
            vdw = None
    else:
        vdw = None

    if comp_parameters.get('functional') in SCAN_list:
        uis['METAGGA'] = 'R2SCAN'
        uis['ALGO'] = 'All'
        uis['LELF'] = False  # otherwise KPAR >1 crashes

    if static_type.endswith('follow_up'):
        uis['ISTART'] = 1
        uis['LREAL'] = '.FALSE.'
        uis['NELMDL'] = -1
    elif static_type.endswith('nsfc'):
        uis['ISTART'] = 1
        uis['LREAL'] = '.FALSE.'
        uis['ICHARG'] = 11
        uis['NELMDL'] = -1

    if 'kspacing' in comp_parameters:
        uis['KSPACING'] = comp_parameters['kspacing']
        uis['KGAMMA'] = True
        kpoints = None
    elif 'k_dens' in comp_parameters:
        if static_type.startswith('slab_') or static_type.startswith('interface_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDensity(structure,
                               comp_parameters['k_dens'],
                               is_slab=is_slab,
                               force_gamma=True)
        kpoints = KPTS.get_kpoints()
    else:
        if static_type.startswith('slab_') or static_type.startswith('interface_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDensity(structure,
                               k_dens_default,
                               is_slab=is_slab,
                               force_gamma=True)
        kpoints = KPTS.get_kpoints()
    uks = kpoints

    if comp_parameters.get('functional') == 'LDA':
        upf = 'LDA_54'
    else:
        upf = 'PBE_54'
        
    if static_type.startswith('bulk_epsilon'):
        uis['LEPSILON'] = True
        uis['KPAR'] = 2
        uis['IBRION'] = 8
        if comp_parameters.get('is_metal', False):
            uis['LPEAD'] = False
        else:
            uis['LPEAD'] = True

    if comp_parameters.get('functional') in SCAN_list:
        vis = MPScanStaticSet(structure, user_incar_settings = uis, vdw = vdw,
                              user_kpoints_settings = uks,
                              user_potcar_functional = 'PBE_54',
                              user_potcar_settings={'W': 'W_sv'})
    else:
        vis = MPStaticSet(structure, user_incar_settings = uis, vdw = vdw,
                          user_kpoints_settings = uks,
                          user_potcar_functional = upf,
                          user_potcar_settings={'W': 'W_sv'})
        
    return vis


def get_custom_vasp_relax_settings(structure, comp_parameters, relax_type,
                                   k_dens_default=8.5, ups={'W': 'W_sv'}):
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
    k_dens_default : float, optional
        Specifies the default kpoint density if no k_dens or kspacing key
        is found in the comp_parameters dictionary. The default is 12.5
    ups : dict, optional
        Specify user potcar settings, e.g. which potcar to choose for which
        element. The default is there to fix an issue with tungsten, where
        MP uses W_pv, which is depreciated by VASP and replaced with W_sv.
        The default is {'W': 'W_sv'}

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
                     'bulk_shape_relax', 'bulk_pos_shape_relax',
                     'slab_shape_relax', 'slab_pos_relax',
                     'interface_shape_relax', 'interface_pos_relax',
                     'interface_z_relax']

    if relax_type not in allowed_types:
        raise SystemExit('relax type is not known. Please select from: {}'
                         .format(allowed_types))

    SCAN_list = ['scan', 'rscan', 'r2scan', 'Scan', 'Rrscan', 'R2scan',
                 'SCAN', 'RSCAN', 'R2SCAN']

    # Set user incar settings:
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
    uis['SYMPREC'] = 1e-04  # compat some issues that VASP 6.2 has with kpoint lattices
    #include electrons with l-quantum number up to 6 into the mixer. Helps with convergence
    #maybe include functionality that sets LMAXMIX dependent on periodic table group
    # e.g.:
    # for element in structure.composition.elements:
    #     if element.group == 3:
    #         uis['LMAXMIX'] = 6
    #     elif element.group in [4,5,6,7,8,9,10,11,12]:
    #         uis['LMAXMIX'] = 4
    uis['LMAXMIX'] = 6 

    if structure.num_sites < 20:
        uis['LREAL'] = '.FALSE.'

    # Adjust mixing for slabs that have a very large c axis:
    if structure.lattice.matrix[-1, 1] > 50.0:
        uis['AMIN'] = 0.05

    if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
        uis['NELMDL'] = -15
        uis['EDIFFG'] = -0.015
        uis['NELM'] = 200
        # Use a slightly slower but more stable algorithm for the electrons
        uis['ALGO'] = 'Normal'
        # Turn on linear mixing
        # uis['AMIX'] = 0.2
        # uis['BMIX'] = 0.0001
        # uis['AMIX_MAG'] = 0.8
        # uis['BMIX_MAG'] = 0.0001
        # set dipole corrections.
        try:
            uis['DIPOL'] = list(structure.center_of_mass)
            uis['IDIPOL'] = 3
            uis['EPSILON'] = comp_parameters.get('epsilon') or 1.0
        except:
            uis['IDIPOL'] = 3
            uis['EPSILON'] = comp_parameters.get('epsilon') or 1.0
    else:
        uis['NELMDL'] = -6
        uis['EDIFFG'] = -0.01
        uis['NELM'] = 100
        uis['ALGO'] = 'Fast'

    if relax_type.startswith('bulk_'):
        uis['IBRION'] = 1

    if relax_type.endswith('full_relax'):
        uis['ISIF'] = 3
    elif relax_type.endswith('pos_relax'):
        uis['ISIF'] = 2
    elif relax_type.endswith('z_relax'):
        uis['ISIF'] = 2
        # Set up selective dynamics array for the structrues site property
        sd_array = []
        for i in range(len(structure.sites)):
            sd_array.append([False, False, True])
        structure.add_site_property('selective_dynamics', sd_array)
    elif relax_type.endswith('vol_relax'):
        uis['ISIF'] = 7
    elif relax_type.endswith('pos_shape_relax'):
        uis['ISIF'] = 4
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
            uis['ISMEAR'] = 1
        elif relax_type.startswith('bulk'):
            uis['SIGMA'] = 0.05
            uis['ISMEAR'] = -5
        else:
            uis['SIGMA'] = 0.1
            uis['ISMEAR'] = 0
    else:
        uis['SIGMA'] = 0.1
        uis['ISMEAR'] = 0

    # set van der Waals functional. Note that as of now, 'functional' must be
    # specified for vdw to work!
    if set(('use_vdw', 'functional')) <= comp_parameters.keys():
        if comp_parameters['use_vdw']:
            if comp_parameters.get('functional') in SCAN_list:
                vdw = 'rVV10'
            else:
                vdw = 'optB86b'
        else:
            vdw = None
    else:
        vdw = None

    if 'kspacing' in comp_parameters:
        uis['KSPACING'] = comp_parameters['kspacing']
        uis['KGAMMA'] = True
        kpoints = None
    elif 'k_dens' in comp_parameters:
        if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDensity(structure,
                               comp_parameters['k_dens'],
                               is_slab=is_slab,
                               force_gamma=True)
        kpoints = KPTS.get_kpoints()
    else:
        if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDensity(structure,
                               k_dens_default,
                               is_slab=is_slab,
                               force_gamma=True)
        kpoints = KPTS.get_kpoints()
    uks = kpoints

    if comp_parameters.get('functional') == 'LDA':
        upf = 'LDA_54'
    else:
        upf = 'PBE_54'

    if 'functional' in comp_parameters:
        if comp_parameters.get('functional') in SCAN_list:
            # Algo All does not play well with tetrahedron method
            if 'is_metal' in comp_parameters:
                if not comp_parameters['is_metal']:
                    uis['SIGMA'] = 0.1
                    uis['ISMEAR'] = 0
            uis['METAGGA'] = 'R2SCAN'
            uis['ALGO'] = 'All'
            uis['LELF'] = False #otherwise KPAR >1 crashes
            vis = MPScanRelaxSet(structure, user_incar_settings = uis,
                                vdw = vdw, user_kpoints_settings = uks,
                                user_potcar_functional = 'PBE_54',
                                user_potcar_settings={'W': 'W_sv'})
        else:
            vis = MPRelaxSet(structure, user_incar_settings = uis, vdw = vdw,
                            user_kpoints_settings = uks,
                            user_potcar_functional = 'PBE_54',
                            user_potcar_settings={'W': 'W_sv'})
    else:
        vis = MPRelaxSet(structure, user_incar_settings = uis, vdw = vdw,
                        user_kpoints_settings = uks,
                        user_potcar_functional = upf,
                        user_potcar_settings={'W': 'W_sv'})

    return vis
