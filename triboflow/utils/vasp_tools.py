import os
import subprocess
import math

from pymatgen.io.vasp.inputs import Kpoints
from kpLib.interface import get_kpoints as get_generalized_kpoints
from pymatgen.io.vasp.sets import (MPRelaxSet, MPScanRelaxSet, MPStaticSet,
                                   MPScanStaticSet)

from triboflow.utils.file_manipulation import remove_matching_files

class MeshFromDenisty:
    """
    Class to find classic or generalized Monkhorst-Pack meshes which may
    or my not be Gamma-Centred from a given k-point density.
    Provides also capabilites to check if meshes generated from similar
    densities are equivalent. Meshes for slabs are adapted according to
    the lower symmetry conditions.
    
    GENERALIZED MONHORST-PACK MESHES NOT YET FULLY TESTED. SHOULD NOT BE USED!
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
        self.klm  = structure.lattice.reciprocal_lattice.abc
        
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
        k1 = math.ceil(k*density)
        k2 = math.ceil(l*density)
        k3 = math.ceil(m*density)       
        if self._is_slab():
            k3 = 1
        return (k1, k2, k3)
    
    def __make_generalized_mesh(self, density):
        """Return a generalized Monkhorst-Pack mesh using kpLib.
        
        More information can be found here:
        http://muellergroup.jhu.edu/K-Points.html
        https://pypi.org/project/kpLib/

        Parameters
        ----------
        density : float
            The minimum allowed distance between lattice points on the
            real-space superlattice (rmin in the paper). This determines the
            density of the k-point grid.

        Returns
        -------
        kpoints : pymatgen.io.vasp.Kpoints
            Generalized Monkhorst-Pack mesh.
        k_dict : dict
            dictionary with information about the Kpoints object.

        """
        
        k_dict = get_generalized_kpoints(structure=self.struct,
                                               include_gamma=None,
                                               symprec=1e-5,
                                               use_scale_factor=False,
                                               minDistance=density,
                                               minTotalKpoints=4)
        
        d_coords = []
        d_weights = []
        
        for i, w in enumerate(k_dict['weights']):
            if w > 0:
                d_coords.append(k_dict['coords'][i])
                d_weights.append(w)
        
        kpoints = Kpoints(comment='GMP-M; '
                          'Rmin = {}; '
                          'TotPoints = {}; '
                          'DistinctPoints = {}'.format(
                              k_dict['min_periodic_distance'],
                              k_dict['num_total_kpts'],
                              k_dict['num_distinct_kpts']),
                          num_kpts=k_dict['num_distinct_kpts'],
                          style=Kpoints.supported_modes.Reciprocal,
                          kpts=d_coords,#[:k_dict['num_distinct_kpts']],
                          kpts_weights=d_weights)
        return kpoints, k_dict
       
    def _is_slab(self):
        """Figures out if the passed structure is a slab.
        
        Automatic detection might fail for slabs that are set up in a non
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
        #has_odd = any(i % 2 == 1 for i in mesh)
        
        if is_hexagonal or self.force_gamma:
            kpoints = Kpoints.gamma_automatic(kpts=mesh)
        else:
            kpoints = Kpoints.monkhorst_automatic(kpts=mesh)
        
        return kpoints
    
    def get_generalised_kpoints(self):
        """Returns a generalized Monkhorst-Pack mesh as a Kpoints object.
        

        Returns
        -------
        kpoints : pymatgen.io.vasp.Kpoints
            Generalized Monkhorst-Pack mesh.

        """
        
        kpoints, kpoints_dict = self.__make_generalized_mesh(self.dens)
        
        return kpoints
    
    def are_generalised_meshes_the_same(self):
        """Compares list of kpoints for generalized k-meshes.
        
        To test if a different target density actually provides a different
        mesh than a reference density.
        

        Returns
        -------
        bool
            True if meshes are the same, False otherwise.

        """
        kpoints_1, info_1 = self.__make_generalized_mesh(self.dens)
        kpoints_2, info_2 = self.__make_generalized_mesh(self.compare_dens)
        
        coords_1 = info_1['coords']
        coords_2 = info_2['coords']
        
        if coords_1 == coords_2:
            return True
        else:
            return False
    
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
        

def get_emin(potcar):
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

    precalc = ['MINDISTANCE = '+str(k_dist),
                'MINTOTALKPOINTS = 4',
                'GAPDISTANCE = 6 ',
                'MONOCLINIC_SEARCH_DEPTH = 2500',
                'TRICLINIC_SEARCH_DEPTH = 1500']
    if Vasp6:
        precalc.append('WRITE_LATTICE_VECTORS = True')
                
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
    remove_matching_files(['KPOINTS*', 'POSCAR*', 'INCAR', 'PRECALC'])

    return KPTS

def get_custom_vasp_static_settings(structure, comp_parameters, static_type):
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
    #uis['SYMPREC'] = 1e-06
    
    if static_type.endswith('from_scratch'):
        uis['ICHARG'] = 2
        uis['LAECHG'] = '.FALSE.'
    
    if structure.num_sites < 20:
        uis['LREAL'] = '.FALSE.'
        
    #Adjust mixing for slabs that have a very large c axis:
    if structure.lattice.matrix[-1,1] > 50.0:
        uis['AMIN'] = 0.05
        
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
            else:
                vdw = 'optB86b'
        else:
            vdw = None
    else:
        vdw = None
    
    if comp_parameters.get('functional') == 'SCAN':
        uis['METAGGA'] = 'SCAN'
        uis['ALGO'] = 'All'
        uis['LELF'] = False #otherwise KPAR >1 crashes
        
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
        #kpoints = Kpoints.automatic(structure, comp_parameters['k_dens'])
        if static_type.startswith('slab_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDenisty(structure,
                               comp_parameters['k_dens'],
                               is_slab=is_slab,
                               force_gamma=True)
        kpoints = KPTS.get_kpoints()
        
    else:
        #if no k-density is supplied in the comp_parameters, use a large value
        #kpoints = Kpoints.automatic(structure, 50)
        if static_type.startswith('slab_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDenisty(structure, 50, is_slab=is_slab, force_gamma=True)
        kpoints = KPTS.get_kpoints()
    
    #If a structure has a vacuum layer, set the third kpoints devision to 1
    #by force.
    # if static_type.startswith('slab_'):
    #     uks = Kpoints(comment=kpoints.comment+'  adjusted for slabs',
    #                 num_kpts=0,
    #                 kpts=[[kpoints.kpts[0][0], kpoints.kpts[0][1], 1]])
    # else:
    #     uks = kpoints
    uks = kpoints
    
    if comp_parameters.get('functional') == 'SCAN':
        vis = MPScanStaticSet(structure, user_incar_settings = uis, vdw = vdw,
                              user_kpoints_settings = uks,
                              user_potcar_functional = 'PBE_54')
    else:
        vis = MPStaticSet(structure, user_incar_settings = uis, vdw = vdw,
                          user_kpoints_settings = uks,
                          user_potcar_functional = 'PBE_54')
        
    return vis
      
def get_custom_vasp_relax_settings(structure, comp_parameters, relax_type):
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
                     'bulk_shape_relax', 'slab_shape_relax', 'slab_pos_relax',
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
    #uis['SYMPREC'] = 1e-06
    
    if structure.num_sites < 20:
        uis['LREAL'] = '.FALSE.'
        
    #Adjust mixing for slabs that have a very large c axis:
    if structure.lattice.matrix[-1,1] > 50.0:
        uis['AMIN'] = 0.05
        
    
    if relax_type.startswith('slab_') or relax_type.startswith('interface_'):
        uis['NELMDL'] = -15
        uis['EDIFFG'] = -0.015
        uis['NELM'] = 200
        #Use a slightly slower but more stable algorithm for the electrons
        uis['ALGO'] = 'Normal'
        # Turn on linear mixing
        # uis['AMIX'] = 0.2
        # uis['BMIX'] = 0.0001
        # uis['AMIX_MAG'] = 0.8
        # uis['BMIX_MAG'] = 0.0001
        
    else:
        uis['NELMDL'] = -6
        uis['EDIFFG'] = -0.01
        uis['NELM'] = 100
    
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
            uis['ISMEAR'] = 1
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
        #kpoints = Kpoints.automatic(structure, comp_parameters['k_dens'])
        if relax_type.startswith('slab_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDenisty(structure,
                               comp_parameters['k_dens'],
                               is_slab=is_slab,
                               force_gamma=True)
        kpoints = KPTS.get_kpoints()
        
    else:
        #if no k-density is supplied in the comp_parameters, use a large value
        #kpoints = Kpoints.automatic(structure, 5000)
        if relax_type.startswith('slab_'):
            is_slab = True
        else:
            is_slab = False
        KPTS = MeshFromDenisty(structure,
                               50,
                               is_slab=is_slab,
                               force_gamma=True)
        kpoints = KPTS.get_kpoints()
    uks = kpoints
        
    if 'functional' in comp_parameters:
        if comp_parameters['functional'] == 'SCAN':
            #Algo All does not play well with tetrahedron method
            if 'is_metal' in comp_parameters:
                if not comp_parameters['is_metal']:
                    uis['SIGMA'] = 0.1
                    uis['ISMEAR'] = 0
            uis['METAGGA'] = 'SCAN'
            uis['ALGO'] = 'All'
            uis['LELF'] = False #otherwise KPAR >1 crashes
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
