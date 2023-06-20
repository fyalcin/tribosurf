import itertools
from collections import defaultdict, Counter
from typing import Union, Any

import numpy as np
from numpy import ndarray
from pymatgen.analysis.local_env import BrunnerNN_real, VoronoiNN
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.core import Structure, Lattice
from pymatgen.core.interface import Interface
from pymatgen.core.surface import center_slab, Slab, get_symmetrically_distinct_miller_indices, SlabGenerator
from pymatgen.io.cif import CifParser
from pymatgen.transformations.standard_transformations import SupercellTransformation
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from surfen.utils.misc_tools import attr_to_dict, get_subset_indices, get_pmg_sg_params
from surfen.utils.misc_tools import check_input


class Shaper:

    @staticmethod
    def get_layer_spacings(struct: Union[Structure, Slab, Interface],
                           tol: float = 0.1, direction: int = 2) -> np.array:
        """
        Simple method to calculate the projected heights of the spacings
        between layers in the given structure.

        Parameters
        ----------
        struct : pymatgen.core.surface.Slab
            Pymatgen object to store slabs. Note that the slab
            should be oriented in such a way that the surface should
            be parallel to the plane that the first 2 lattice vectors lie on.
        tol : float, optional
            Tolerance parameter to cluster sites into layers. The default is 0.1.
        direction : int, optional
            Direction in which we calculate the projected height. Allowed
            values are 0, 1, and 2, which correspond to the first, second,
            and third lattice vectors respectively.
            The default is 2.

        Returns
        -------
        list
            list of floats representing the projected distances between layers
            along the surface normal direction in angstroms

        """

        # Layer info that contains the c-coordinates and sites
        layers = Shaper.get_layers(struct, tol, direction)

        # Only the c-coordinates of the layers are needed
        layers_c = sorted(layers.keys())

        # Spacing between consecutive layers are calculated
        d = [x - layers_c[i - 1] for i, x in enumerate(layers_c)]

        # When a periodic boundary is passed, layers wrap over and we get a
        # negative spacing, to correct, we add 1 to negative spacing values
        d = [s + int(s < 0) for s in d]

        # For slabs with the third lattice vector not along miller
        # direction, we need the projected height to also project the vacuum
        # height
        proj_height = Shaper.get_proj_height(struct)

        return np.round([spacing * proj_height for spacing in d], 10)

    @staticmethod
    def get_proj_height(struct: Union[Structure, Slab, Interface], region: str = 'cell', min_vac: float = 4.0,
                        direction: int = 2, tol: float = 0.1) -> float:
        """
        Internal method to calculate the projected height of a specific region.
        For more than one slab region, the total height is calculated.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Main object in pymatgen to store structures.
        region : str, optional
            Region to calculate the projected height for. Can take values
            'cell', 'vacuum', or 'slab'. The default is 'cell'.
        min_vac : float, optional
            Thickness threshold in angstroms to define a region as a
            vacuum region. The default is 4.0.
        direction : int, optional
            Direction in which we calculate the projected height. Allowed
            values are 0, 1, and 2, which correspond to the first, second,
            and third lattice vectors respectively.
            The default is 2.
        tol : float, optional
            Tolerance parameter to cluster sites into layers. The default is 0.1.

        Raises
        ------
        ValueError
            Simple check for region keyword to see if it's one of allowed values.

        Returns
        -------
        proj_height : float
            Projected height of the region. The thickness or height is projected
            along the hkl direction which is assumed to be the normal direction
            to the first two lattice vectors of the passed structure.

        """

        # proj_height = Shaper.get_hkl_projection(struct.lattice.matrix[direction], struct)
        lateral_vecs = [struct.lattice.matrix[i] for i in range(3) if i != direction]
        normal = np.cross(*lateral_vecs)
        normal /= np.linalg.norm(normal)
        vec_to_proj = struct.lattice.matrix[direction]
        proj_height = np.abs(np.dot(vec_to_proj, normal))
        if region == "cell":
            return proj_height
        elif region == "slab" or region == "vacuum":
            spacings = Shaper.get_layer_spacings(struct, tol)
            slab_height = sum([s for s in spacings if s < min_vac])
            return slab_height if region == "slab" else proj_height - slab_height
        else:
            raise ValueError('Region must be one of "cell", "vacuum", or "slab"')

    @staticmethod
    def resize(struct: Union[Structure, Slab, Interface], slab_thickness: Union[int, float] = None,
               vacuum_thickness: Union[int, float] = None, tol: float = 0.1,
               chunk_size: int = 1, min_thick_A: Union[int, float] = None, center: bool = True,
               min_vac: Union[int, float] = 4.0) -> Union[Structure, Slab, Interface]:
        """
        resize the input slab with the desired slab thickness in
        number of layers and the vacuum region in Angstroms. All the attributes
        of sites are preserved by the resizing process.

        Parameters
        ----------
        min_vac : float, optional
            Thickness threshold in angstroms to define a region as a
            vacuum region. The default is 4.0 A.
        struct : pymatgen.core.structure.Structure
            Structure object that is to be resized. Input object
            is not modified with this method.
        slab_thickness : int
            Desired slab thickness in number of layers. Layers will
            be removed from the bottom until the desired thickness is
            reached.
        vacuum_thickness : float
            Desired vacuum region thickness in Angstroms. Lattice
            parameters are modified in order to get the correct vacuum.
        tol : float, optional
            Tolerance value for layering of sites in Angstroms. Used in counting
            of layers. The default is 0.1.
        chunk_size : int, optional
            Number of layers that are removed at once. Used to preserve terminations.
            The default is 1.
        min_thick_A : float, optional
            Minimum slab thickness in Angstroms. If a float is passed, the number of
            layers to remove will be adjusted so that the resized slab is at least
            min_thick_A thick in Angstroms.
            The default is None.
        center : bool, optional
            Whether to center the resized slab between the vacuum region.
            The default is True.

        Returns
        -------
        resized_struct : pymatgen.core.structure.Structure
            resized structure with the desired parameters.

        """
        # Input slab is first centered for the cases where the slab spills
        # outside the box from the top and the bottom
        struct_centered = center_slab(struct.copy())
        initial_thickness = Shaper.get_proj_height(struct_centered, 'slab', min_vac)

        if slab_thickness:
            # Layers (containing sites) are removed from the bottom until
            # the desired slab_thickness is reached
            num_layers = len(Shaper.get_layers(struct_centered, tol))
            layers_to_remove = int(chunk_size * np.floor((num_layers - slab_thickness) / chunk_size))
            if min_thick_A:
                spacings = [spacing for spacing in Shaper.get_layer_spacings(struct_centered, tol) if spacing < 4.0]
                if initial_thickness > min_thick_A:
                    while initial_thickness - sum(spacings[:layers_to_remove]) < min_thick_A:
                        layers_to_remove -= chunk_size
                else:
                    print(f'Slab is already smaller than min_thick_A, resizing halted..')
                    layers_to_remove = 0
            if layers_to_remove > 0:
                struct_resized = Shaper.remove_layers(struct_centered, layers_to_remove,
                                                      tol=tol, method='layers')
            else:
                struct_resized = struct_centered
        else:
            struct_resized = struct_centered
        # Vacuum region is modified to the desired thickness
        if vacuum_thickness:
            resized_struct = Shaper.modify_vacuum(struct_resized, vacuum_thickness, min_vac=min_vac)
        else:
            resized_struct = struct_resized

        if not slab_thickness and not vacuum_thickness:
            print(f'Warning! You chose to keep the slab and vacuum thicknesses as they are'
                  'during resize. Make sure this is what you want.')
            resized_struct = struct_centered if center else struct

        # bbs = kwargs.get('bbs')
        # if bbs:
        #     layers_initial = len(Shaper.get_layers(struct_centered, tol))
        #     layers_resized = len(Shaper.get_layers(resized_struct, tol))
        #     diff = layers_initial - layers_resized
        #     shifts = list(bbs.keys())
        #     top_shift = np.round(resized_struct.shift, 4)
        #     top_shift_index = shifts.index(top_shift)
        #     bot_shift = shifts[(top_shift_index - diff) % len(shifts)]
        #     top_bvs = bbs[top_shift]
        #     bot_bvs = bbs[bot_shift]
        #     resized_struct.energy = {'top': top_bvs, 'bottom': bot_bvs}
        return resized_struct

    @staticmethod
    def modify_vacuum(struct: Union[Structure, Slab, Interface], vac_thick: Union[int, float], method: str = 'to_value',
                      center: bool = True, min_vac: Union[int, float] = 4.0) -> Union[Structure, Slab, Interface]:
        """
        Method to modify the vacuum region in a structure.

        Parameters
        ----------
        min_vac : float, optional
            Thickness threshold in angstroms to define a region as a
            vacuum region. The default is 4.0 A.
        struct : pymatgen.core.structure.Structure
            Main object in pymatgen to store structures.
        vac_thick : float
            Vacuum adjustment amount in Angstroms.
        method : str, optional
            Whether to set the vacuum to the desired value or adjust the
            vacuum in the structure by the given value.
            The default is 'to_value'.
        center : bool, optional
            Whether to center the slab in the resulting structure inside
            the vacuum region.
            The default is True.

        Returns
        -------
        modified_struct : pymatgen.core.structure.Structure
            Modified pymatgen Structure object.

        """

        # Check if a Slab or Structure is passed and proceed accordingly
        if 'miller_index' in vars(struct):
            # Necessary slab attributes to resize the Slab
            attrs = ["species", "miller_index", "oriented_unit_cell",
                     "shift", "scale_factor", "reorient_lattice",
                     "reconstruction", "site_properties", "energy"]
            struct_params = attr_to_dict(struct, attrs)
            out_object = Slab
        else:
            # Necessary structure attributes to resize the Structure
            attrs = ["species", "site_properties"]
            struct_params = attr_to_dict(struct, attrs)
            out_object = Structure

        # To avoid issues with fractional coordinates when scaling vacuum,
        # cartesian coordinates are used
        corrected_params = {'coords': struct.cart_coords,
                            'coords_are_cartesian': True}
        struct_params.update(corrected_params)

        # Lattice parameters are generated in order to be modified
        lat_attrs = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        lat_params = attr_to_dict(struct.lattice, lat_attrs)

        # latvec = struct.lattice.matrix
        proj_height = Shaper.get_proj_height(struct)

        # 'c' parameter of the Lattice is modified to adjust vacuum
        # to the desired thickness
        if method == 'to_value':
            initial_vac = Shaper.get_proj_height(struct, 'vacuum', min_vac)
            lat_params['c'] += (vac_thick - initial_vac) * lat_params['c'] / proj_height
        elif method == 'by_value':
            lat_params['c'] += vac_thick * lat_params['c'] / proj_height

        new_lat = Lattice.from_parameters(**lat_params)

        modified_struct = center_slab(out_object(new_lat, **struct_params)) \
            if center else out_object(new_lat, **struct_params)

        return modified_struct

    # @staticmethod
    # def get_hkl_projection(vector, struct):
    #     """
    #     Simple method to calculate the norm of the projection of a vector
    #     along the hkl direction which is assumed to be normal to the plane
    #     formed by the first two lattice vectors of the passed structure.
    #     Useful for structures where the third lattice vector is not in
    #     the same direction as the surface normal.
    #
    #     Parameters
    #     ----------
    #     struct : pymatgen.core.structure.Structure
    #         Main object in pymatgen to store structures.
    #
    #     Returns
    #     -------
    #     float
    #         Projected height of the given structure in the direction
    #         that is normal to the x-y plane
    #
    #     """
    #     latvec = struct.lattice.matrix
    #     normal = np.cross(latvec[0], latvec[1])
    #     normal /= np.linalg.norm(normal)
    #     return np.abs(np.dot(vector, normal))

    @staticmethod
    def get_layers(struct: Union[Structure, Slab, Interface], tol: float = 0.1,
                   direction: int = 2) -> dict:
        """
        Finds the layers in the structure taking z-direction as the primary
        direction such that the layers form planes parallel to xy-plane.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Main object in pymatgen to store structures. Has to be given in a
            way that the first two lattice vectors lie on a plane perpendicular
            to a given miller direction.

        tol : float, optional
            Tolerance parameter to cluster sites into layers. The default is 0.1.

        direction : int, optional
            Which direction to count the layers in. Allowed values are 0, 1, and 2,
            and they correspond to the first, second, and the third lattice
            vectors respectively.
            The default is 2.

        Returns
        -------
        layers : dict
            Dictionary with keys as z-coords of layers and values as the
            indices of sites that belong to that layer.

        """
        # number of sites in the structure
        n = len(struct)
        frac_coords = struct.frac_coords

        # initiate a num_sites dimensional square distance matrix and populate
        dist_matrix = np.zeros((n, n))
        for i, j in itertools.combinations(list(range(n)), 2):
            if i != j:
                cdist = frac_coords[i][direction] - frac_coords[j][direction]
                # cdist = abs(cdist - round(cdist)) * proj_height
                cdist = abs(cdist - round(cdist)) * struct.lattice.abc[direction]
                dist_matrix[i, j] = cdist
                dist_matrix[j, i] = cdist

        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)

        # cluster the sites in the structure based on their c-coordinate
        # and a given tolerance
        clusters = fcluster(z, tol, criterion="distance")
        layers = defaultdict(list)
        for i, v in enumerate(clusters):
            layers[v].append(i)

        # for each layer, find sites that belong to it and assign the first
        # site's c-coord as the c-coord of the layer
        layers = {struct.frac_coords[v[0]][direction]: v for k, v in layers.items()}
        return layers

    @staticmethod
    def remove_layers(slab: Union[Slab, Interface], num_layers: int, tol: float = 0.1, method: str = 'target',
                      position: str = 'bottom', center: bool = True) -> Union[Slab, Interface]:
        """
        Removes layers from the bottom of the slab while updating the number
        of bonds broken in the meantime.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Standard pymatgen Slab object.
        num_layers : int
            Number of layers to remove from the structure
        tol : float, optional
            Tolerance to use in the identification of the layers.
            The default is 0.1.
        method : str, optional
            Whether to remove num_layers or remove layers until the
            structure has num_layers number of layers in total.
            Options are 'target' and 'layers'. The default is 'target'.
        position : string, optional
            Side on which the sites should be removed.
            Available options are 'top' and 'bottom'. The default is 'bottom'.
        center : bool, optional
            Whether to center the slab in the vacuum after removing layers.
            The default is 'True'.

        Returns
        -------
        slab_copy : pymatgen.core.surface.Slab
            Copy of the input Slab structure with layers removed.

        """
        layers = Shaper.get_layers(slab, tol)
        if num_layers > len(layers):
            raise ValueError('Number of layers to remove/target can\'t exceed \
                             the number of layers in the given slab.')
        c_coords = sorted(layers.keys())
        if method == "layers":
            to_remove = c_coords[:num_layers] if position == "bottom" \
                else c_coords[len(c_coords) - num_layers:]
        elif method == "target":
            to_remove = c_coords[:len(c_coords) - num_layers] if position == "bottom" \
                else c_coords[num_layers:]
        else:
            raise ValueError(f'{method} is not a valid method. Please use either "layers" or "target".')
        indices_list = [layers[c_coord] for c_coord in to_remove]
        flat_list = [item for sublist in indices_list for item in sublist]
        slab_copy = slab.copy()
        slab_copy.remove_sites(flat_list)
        return center_slab(slab_copy) if center else slab_copy

    @staticmethod
    def get_average_layer_spacing(slab: Union[Slab, Interface], tol: float = 0.1,
                                  vacuum_threshold: Union[int, float] = 6.0) -> ndarray:
        """
        Compute the average distance between the slabs layers disregarding the
        vacuum region.

        Parameters
        ----------
        tol : float, optional
            Tolerance to use in the identification of the layers.
            The default is 0.1 Angstroms.
        slab : pymatgen.core.surface.Slab
            Standard pymatgen Slab object.
        vacuum_threshold : float, optional
            Regions larger than this will be treated as vacuum and will not be
            treated as an interlayer spacing. The default is 6.0
        Returns
        -------
        av_spacing : float
            Average layer spacing

        """
        spacings = Shaper.get_layer_spacings(slab, tol)
        spacings_no_vac = np.delete(spacings,
                                    np.where(spacings >= vacuum_threshold))
        av_spacing = np.mean(spacings_no_vac)
        return av_spacing

    @staticmethod
    def get_bonds(struct: Union[Slab, Structure], method: str = 'covalent_radii',
                  dtol: float = 0.20) -> dict:
        """
        Finds all unique bonds in the structure and orders them by bond strength
        using bond valance method and with the assumption that the ideal bond length
        = CovalentRadius(site1) + CovalentRadius(site2)

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Conventional standard structure that is used to generate the slabs.

        method : string, optional
            Method used to calculate the bond valence parameters
            - 'covalent_radii' : sets the 'ideal' bond length R_0 as the sum of
                the covalent radii of bonded atoms
            - 'BVparams' : uses a .cif file of a list of fitted bond valence
                parameters from various sources.

        dtol : float, optional
            Added tolerance to form a bond for the dictionary passed to the
            slab generation algorithm.

        Returns
        -------
        dict : Collection of bonds that has a 'weight' within a delta of the highest
        weight.
        """
        # struct = struct.get_primitive_structure()
        bnn = BrunnerNN_real(cutoff=2 * max(struct.lattice.abc))
        # cutoff = 1.5*max(struct.lattice.abc)
        species, indices = np.unique([str(x) for x in struct.species],
                                     return_index=True)
        bonds = {}
        wmax = 0
        for i, site_index in enumerate(indices):
            sp1 = species[i]
            # for neighbor in struct.get_neighbors(site=struct[site_index], r=cutoff):
            for neighbor in bnn.get_nn_info(struct, site_index):
                neighbor = neighbor['site']
                sp2 = str(neighbor.specie)
                dist = neighbor.nn_distance
                if method == 'covalent_radii':
                    cr = CovalentRadius().radius
                    R_0 = cr[sp1] + cr[sp2]
                    b = 0.37
                elif method == 'BVparams':
                    ciffile = CifParser('bvparm2020.cif').as_dict()
                    a1 = ciffile['BOND_VALENCE_PARAMETERS_2020-11-25']
                    b1 = a1['_valence_param_atom_1']
                    b2 = a1['_valence_param_atom_2']
                    b3 = a1['_valence_param_Ro']
                    b4 = a1['_valence_param_B']
                    c = defaultdict(list)
                    for i in range(len(b1)):
                        try:
                            c[(b1[i], b2[i])].append((b3[i], b4[i]))
                        except:
                            c[(b1[i], b2[i])] = [(b3[i], b4[i])]
                    # R_0 = sum([float(i[0]) for i in c[(sp1,sp2)]])/len(c[(sp1,sp2)])
                    try:
                        R_0 = float(c[(sp1, sp2)][0][0])
                        b = float(c[(sp1, sp2)][0][1])
                    except:
                        cr = CovalentRadius().radius
                        R_0 = cr[sp1] + cr[sp2]
                        b = 0.37
                elif not isinstance(method, str):
                    raise TypeError("method argument must be a string")
                else:
                    raise ValueError('method can either be covalent_radii or BVparams')
                w = np.exp((R_0 - dist) / b)
                wmax = w if w > wmax else wmax
                if ((sp1, sp2) not in bonds) and ((sp2, sp1) not in bonds):
                    bonds[(sp1, sp2)] = (dist * (1 + dtol), w)
        # bonds = {k: v[0] for k, v in bonds.items() if abs(v[1]-wmax)/wmax <= wtol}
        return bonds

    @staticmethod
    def get_c_ranges_old(struct, nn_method='all', cutoff=5.0, weights='bvs', c_range_pmg=True):
        """
        Calculates all the bond valence sums for bonds that would be broken for each
        possible shift in the given structure.

        Parameters
        ----------
        weights : string, optional
            Method used to calculate the bond valence parameters
            - 'bvs' : uses the bond valence sums from the pymatgen BrunnerNN_real
        struct : pymatgen.core.structure.Structure
            Pymatgen Structure object.
        nn_method: str, optional
            Nearest-neighbor algorithm to be used. Currently, supports
            'all' and 'BNN'. For more info, check out BrunnerNN documentation.
        cutoff : float, optional
            Cutoff radius used in neighbor searching. The value is in Angstroms.
            The default value is 5.0A.

        Returns
        -------
        c_ranges : list
            List with elements describing every bond in the structure, with
            the endpoints of the bond and the bond valence values for each bond.

        """
        cr = CovalentRadius().radius
        if nn_method == 'all':
            nn_list = struct.get_all_neighbors(r=cutoff)
        elif nn_method == 'BNN':
            bnn = BrunnerNN_real(cutoff=cutoff)
            nn_list = bnn.get_all_nn_info(struct)
        else:
            raise ValueError('"nn_method" must be one of "all" or "BNN"')

        print('Neighbor finding done!')
        if nn_method == 'all' and weights == 'equal':
            print(f'WARNING! Setting nn_method to "all" and weights to "equal" will lead to incorrect results.')

        c_ranges = []
        species = [str(i) for i in struct.species]
        c_ranges_uc = []
        for s_index, site in enumerate(struct):
            for nn in nn_list[s_index]:
                if nn_method == 'BNN':
                    image = nn['image']
                    nn = nn['site']
                c_range = np.round(sorted([site.frac_coords[2], nn.frac_coords[2]]), 3)
                # c_range = np.round([site.frac_coords[2], nn.frac_coords[2]], 3)

                # if c_range[0] != c_range[1]:
                if True:
                    # if image == (0, 0, 0):
                    #     continue
                    # if c_range.tolist() in c_ranges_uc:
                    #     continue
                    # else:
                    #     c_ranges_uc.append(c_range.tolist())
                    nn_site_index = nn.index
                    sp1 = species[s_index]
                    sp2 = species[nn_site_index]
                    dist = nn.nn_distance
                    if weights == 'bvs':
                        w = Shaper.calculate_bv(cr[sp1], cr[sp2], dist)
                    elif weights == 'equal':
                        w = 1
                    elif weights == 'bvs_squared':
                        w = Shaper.calculate_bv(cr[sp1], cr[sp2], dist) ** 2
                    else:
                        raise Exception(f'Weights for c_ranges must be one of "bvs", "equal", or "bvs_squared".')
                    bv = ((sp1, s_index), (sp2, nn_site_index), dist, w)
                    if c_range_pmg:
                        # while c_range[1] <= 0:
                        #     c_range += 1
                        c_ranges.append((c_range, bv))

                        # if c_range[0] < 0:
                        #     c_ranges.append((0, c_range[1], bv))
                        #     c_ranges.append((c_range[0] + 1, 1, bv))
                        # elif c_range[1] > 1:
                        #     c_ranges.append((c_range[0], 1, bv))
                        #     c_ranges.append((0, c_range[1] - 1, bv))
                        # else:
                        #     c_ranges.append((c_range[0], c_range[1], bv))
                    else:
                        c_ranges.append((c_range[0], c_range[1], bv))
        return c_ranges

    @staticmethod
    def get_c_ranges(struct: Union[Slab, Structure], nn_method: str = 'all', cutoff: Union[int, float] = 5.0):
        """
        Calculates all the bond valence sums for bonds that would be broken for each
        possible shift in the given structure.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Pymatgen Structure object.
        nn_method: str, optional
            Nearest-neighbor algorithm to be used. Currently supports
            'all' and 'BNN'. For more info, check out BrunnerNN documentation.
        cutoff : float, optional
            Cutoff radius used in neighbor searching. The value is in Angstroms.
            The default value is 5.0A.

        Returns
        -------
        c_ranges : list
            List with elements describing every bond in the structure, with
            the endpoints of the bond and the bond valence values for each bond.

        """
        cr = CovalentRadius().radius
        # cutoff = 2 * max([a[0] for a in list(Shaper.get_bonds(struct).values())])
        if nn_method == 'all':
            # print('Finding all neighbors!')
            nn_list = struct.get_all_neighbors(r=cutoff)
        elif nn_method == 'BNN':
            # print('Using BrunnerNN algorithm!')
            bnn = BrunnerNN_real(cutoff=cutoff)
            nn_list = bnn.get_all_nn_info(struct)
        elif nn_method == 'VNN':
            vnn = VoronoiNN(cutoff=cutoff)
            nn_list = vnn.get_all_nn_info(struct)
        else:
            raise ValueError('"nn_method" must be one of "all" or "BNN"')

        # if nn_method == 'all' and weight == 'equal':
        #     print(f'WARNING! Setting nn_method to "all" and weights to "equal" will lead to incorrect results.')

        c_ranges = []
        species = [str(i) for i in struct.species]
        for s_index, site in enumerate(struct):
            for nn in nn_list[s_index]:
                if nn_method == 'BNN':
                    # w = nn['weight']
                    # if np.round(w, 2) < 1.0:
                    #     continue
                    nn = nn['site']
                c_range = np.round(sorted([site.frac_coords[2], nn.frac_coords[2]]), 4)
                sp1, sp2, dist = species[s_index], species[nn.index], nn.nn_distance
                cr_dist = (cr[sp1], cr[sp2], dist)
                c_ranges.append((c_range, cr_dist))
        return c_ranges, nn_list

    @staticmethod
    def calculate_bv(r1: float, r2: float, bond_dist: float) -> float:
        b, r_0 = 0.37, r1 + r2
        return np.exp((r_0 - bond_dist) / b)

    @staticmethod
    def get_surface_area(struct: Union[Structure, Slab, Interface]) -> float:
        mat = struct.lattice.matrix
        return np.linalg.norm(np.cross(mat[0], mat[1]))

    @staticmethod
    def bonds_by_shift_old(sg, nn_methods=('all',), tol=0.1, cutoff=5.0, weights=('bvs',)):
        """
        Calculates the bond valence sums of the broken bonds corresponding to
        all the possible shifts

        Parameters
        ----------
        weights : list, optional
            List of weights to be used for each shift. The default is ['bvs'].
        nn_methods : list, optional
            List of nearest-neighbor algorithms to be used. Currently supports
            'all' and 'BNN'. For more info, check out BrunnerNN documentation.
        sg : pymatgen.core.surface.SlabGenerator
            Pymatgen SlabGenerator object to extract the possible shifts in a
            specific orientation.
        tol : float, optional
            Tolerance value used in the layering of sites in units of Angstroms.
            The default is 0.1.
        cutoff : float, optional
            Cutoff radius used in neighbor searching. The value is in Angstroms.
            The default value is 5.0A.

        Returns
        -------
        bbs : dict
            Dictionary with keys as shifts and values as the bond valence sum of the
            broken bonds at each shift, scaled by the area of the x-y plane.

        """
        ouc = sg.oriented_unit_cell
        shifts = np.round(sg._calculate_possible_shifts(tol=tol), 4)

        scaling_matrix = ((1, 0, 0), (0, 1, 0), (0, 0, 3))
        ouc.make_supercell(scaling_matrix=scaling_matrix)

        shifts = [0.33 + shift / 3 for shift in shifts]

        bbs = {}
        for nn_method in nn_methods:
            c_ranges, nn_list = Shaper.get_c_ranges(ouc, nn_method, cutoff)
            for weight in weights:
                bbs[(nn_method, weight)] = {}
                for shift in shifts:
                    tmp = sum(
                        [Shaper.calculate_bv(*c[1]) if weight == 'bvs' else 1 for c in c_ranges if
                         c[0][0] < shift < c[0][1]])
                    bbs[(nn_method, weight)][shift] = np.round(tmp, 4)
        return bbs

    @staticmethod
    def bonds_by_shift(sg: SlabGenerator, nn_methods: tuple = ('all',), tol: float = 0.1,
                       cutoff: Union[int, float] = 5.0,
                       weights: tuple = ('bvs',)) -> tuple[dict[tuple[tuple, tuple], dict[Any, Any]], dict]:
        """
        Calculates the bond valence sums of the broken bonds corresponding to
        all the possible shifts

        Parameters
        ----------
        nn_methods : list, optional
            List of nearest-neighbor algorithms to be used. Currently supports
            'all' and 'BNN'. For more info, check out BrunnerNN documentation.
        sg : pymatgen.core.surface.SlabGenerator
            Pymatgen SlabGenerator object to extract the possible shifts in a
            specific orientation.
        nn_method: str, optional
            Nearest-neighbor algorithm to be used. Currently supported algorithms are
            'all' and 'BNN'. For more info, check out BrunnerNN documentation.
        tol : float, optional
            Tolerance value used in the layering of sites in units of Angstroms.
            The default is 0.1.
        cutoff : float, optional
            Cutoff radius used in neighbor searching. The value is in Angstroms.
            The default value is 5.0A.
        weights : string, optional
            Weight scheme used when counting the bonds broken. Supported options are
            "bvs", "bvs_squared", and "equal". "bvs" and "bvs_squared" schemes calculate
            and add up the bond valence sum and square of bond valence sum of broken bonds,
            respectively, while "equal" treats each broken bond as equal. Usage of "equal"
            should be avoided when nn_method is set to "all" as this will consider each broken
            bond, even when there is no actual bond due to a large bond length equally.

        Returns
        -------
        bbs : dict
            Dictionary with keys as shifts and values as the bond valence sum of the
            broken bonds at each shift, scaled by the area of the x-y plane.

        """
        if not isinstance(nn_methods, list):
            nn_methods = [nn_methods]
        if not isinstance(weights, list):
            weights = [weights]

        ouc = sg.oriented_unit_cell
        bvs_bulk = {}
        for nn_method in nn_methods:
            cutoff_nn = 2 * cutoff if nn_method == 'BNN' else cutoff
            c_ranges_bulk, nn_list_bulk = Shaper.get_c_ranges(ouc, nn_method, cutoff_nn)
            for weight in weights:
                # if weight == 'bvs' or nn_method == 'BNN':
                #     cutoff_nn = 2 * cutoff
                # else:
                #     cutoff_nn = cutoff
                c_ranges_bulk, nn_list_bulk = Shaper.get_c_ranges(ouc, nn_method, cutoff_nn)
                bv_bulk = sum([Shaper.calculate_bv(*c[1]) if weight == 'bvs' else 1 for c in c_ranges_bulk])
                bv_bulk_per_atom = bv_bulk / ouc.num_sites
                bvs_bulk[(nn_method, weight)] = bv_bulk_per_atom

        slabs = {np.round(slab.shift, 4): slab for slab in sg.get_slabs(ftol=tol)}
        # print(f'slab area for {sg.miller_index} is {[slab.surface_area for slab in slabs.values()]}')

        # for shift, slab in slabs.items():
        #     slab.to('poscar', f'{slab.composition.reduced_formula}_{slab.miller_index}.vasp')

        bbs = {}
        for nn_method in nn_methods:
            cutoff_nn = 2 * cutoff if nn_method == 'BNN' else cutoff
            for weight in weights:
                # if weight == 'bvs' or nn_method == 'BNN':
                #     cutoff_nn = 2 * cutoff
                # else:
                #     cutoff_nn = cutoff
                bbs[(nn_method, weight)] = {}
                for shift, slab in slabs.items():
                    c_ranges_slab, nn_list_slab = Shaper.get_c_ranges(slab, nn_method, cutoff_nn)
                    bv_slab = sum([Shaper.calculate_bv(*c[1]) if weight == 'bvs' else 1 for c in c_ranges_slab])
                    diff = slab.num_sites * bvs_bulk[(nn_method, weight)] - bv_slab
                    layers = Shaper.get_layers(slab)
                    num_sites_top_layer = len(layers[max(layers, key=layers.get)])
                    # bbs[(nn_method, weight)][shift] = np.round((diff / slab.surface_area), 4)
                    bbs[(nn_method, weight)][shift] = (np.round(diff, 4), num_sites_top_layer, slab.surface_area)

        return bbs, slabs

    @staticmethod
    def bonds_by_shift_new(sg: SlabGenerator, bulk_conv: Structure, nn_method: str = 'all',
                           tol: float = 0.1, edge_tol: int = 3, cutoff: Union[int, float] = 5.0
                           ) -> tuple[dict[str, dict[Any, Any]], Any]:
        """
        Calculates the bond valence sums of the broken bonds corresponding to
        all the possible shifts

        Parameters
        ----------
        edge_tol : int, optional
            Number of layers to be considered as edge layers. The default is 3.
        nn_method: str, optional
            Nearest-neighbor algorithm to be used. Currently supported algorithms are
            'all' and 'BNN'. For more info, check out BrunnerNN documentation.
        bulk_conv : pymatgen.core.structure.Structure
            Pymatgen Structure object of the bulk structure.
        sg : pymatgen.core.surface.SlabGenerator
            Pymatgen SlabGenerator object to extract the possible shifts in a
            specific orientation.
        tol : float, optional
            Tolerance value used in the layering of sites in units of Angstroms.
            The default is 0.1.
        cutoff : float, optional
            Cutoff radius used in neighbor searching. The value is in Angstroms.
            The default value is 5.0A.
        weights : string, optional
            Weight scheme used when counting the bonds broken. Supported options are
            "bvs", "bvs_squared", and "equal". "bvs" and "bvs_squared" schemes calculate
            and add up the bond valence sum and square of bond valence sum of broken bonds,
            respectively, while "equal" treats each broken bond as equal. Usage of "equal"
            should be avoided when nn_method is set to "all" as this will consider each broken
            bond, even when there is no actual bond due to a large bond length equally.

        Returns
        -------
        bbs : dict
            Dictionary with keys as shifts and values as the bond valence sum of the
            broken bonds at each shift, scaled by the area of the x-y plane.

        """
        unique_sites = np.unique(bulk_conv.site_properties['bulk_equivalent'], return_index=True)[1]
        c_bulk = {}
        bonds_bulk = {}
        for site_index in unique_sites:
            site = bulk_conv[site_index]
            specie = str(site.specie)
            nn_list = Shaper.get_neighbors(bulk_conv, cutoff, nn_method, site_index)
            c_bulk[site_index] = len(nn_list)
            bonds_bulk[f'{specie}-{site_index}'] = [
                (f'{str(nn.specie)}-{nn.properties["bulk_equivalent"]}',
                 np.round(nn.nn_distance, 6))
                for nn in nn_list]

        slabs = {np.round(slab.shift, 4): slab for slab in sg.get_slabs(ftol=tol)}
        area = list(slabs.values())[0].surface_area

        # for shift, slab in slabs.items():
        #     slab.add_oxidation_state_by_guess()
        #     print(f'Polarity of slab with shift {shift} is {slab.is_polar()}')

        # layer_spacings = [spc for spc in Shaper.get_layer_spacings(list(slabs.values())[0], tol) if spc > 4]
        # print(f'for hkl {sg.miller_index}, {edge_tol} layers is {sum(layer_spacings[-edge_tol:])} thick')

        bbs = {}
        for shift, slab in slabs.items():
            layers = Shaper.get_layers(slab, tol)
            layers_c_sorted = sorted(layers.keys())
            c_top_layer, c_bot_layer = max(layers), min(layers)

            if edge_tol > len(layers_c_sorted) or edge_tol == 99:
                top_layers_c_coords = layers_c_sorted
            else:
                top_layers_c_coords = layers_c_sorted[-edge_tol:]
            # top_layers_c_coords = [a for a in layers_c_sorted if a > c_top_layer - edge_tol]
            # print(f'top layer c coords are {top_layers_c_coords}')

            top_layer_sites = [site for k in top_layers_c_coords for site in layers[k]]
            # top_layer_sites = [item for sublist in top_layer_sites for item in sublist]

            # bot_layers_c_coords = layers_c_sorted[:edge_tol]
            # bot_layers_c_coords = [a for a in layers_c_sorted if a < c_bot_layer + edge_tol]
            # bot_layer_sites = [layers[k] for k in bot_layers_c_coords]
            # bot_layer_sites = [item for sublist in bot_layer_sites for item in sublist]
            # num_sites_top_layer = len(layers[max(layers)])
            # num_sites_bot_layer = len(layers[min(layers)])

            nn_dict = defaultdict(list)
            # c_slab_tot = 0
            # c_bulk_tot = 0
            # n_top = 0
            for site_index in top_layer_sites:
                site = slab[site_index]
                specie = str(site.specie)
                bulk_eq = site.properties['bulk_equivalent']
                nn_list = Shaper.get_neighbors(slab, cutoff, nn_method, site_index)
                nn_list_hr = [(f'{str(nn.specie)}-{nn.properties["bulk_equivalent"]}',
                               np.round(nn.nn_distance, 6)) for nn in nn_list]
                nn_list_hr_bulk = bonds_bulk[f'{specie}-{bulk_eq}']
                diff = Counter(nn_list_hr_bulk) - Counter(nn_list_hr)
                bonds_broken = list(diff.elements())
                if bonds_broken:
                    # n_top += 1
                    # c_slab_tot += c_bulk[nn_method][bulk_eq] - len(bonds_broken)
                    # c_bulk_tot += c_bulk[nn_method][bulk_eq]
                    nn_dict[f'{specie}-{bulk_eq}'] += bonds_broken

            # nn_dict = [aa[-1] if weights[0] == 'bvs' else 1 for sublist in nn_dict.values() for aa in
            #                sublist]

            # bbs[nn_method][shift] = (sum(nn_dict), 1, slab.surface_area)

            bbs[shift] = nn_dict

            # bbs[(nn_method, weight)][shift] = (fact, n_top, slab.surface_area)
            # nn_dict_bot = {}
            # for site_index in bot_layer_sites:
            #     site = slab[site_index]
            #     bulk_eq = site.properties['bulk_equivalent']
            #     nn_list_hr_bulk = bvs_bulk[(nn_method, weight)][bulk_eq]
            #     # nn_list = slab.get_neighbors(site=site, r=cutoff)
            #     nn_list = BNN.get_nn_info(slab, site_index)
            #     nn_list = [nn['site'] for nn in nn_list]
            #     nn_list_hr = [(str(nn.specie), np.round(nn.nn_distance, 8)) for nn in nn_list]
            #     c1 = Counter(nn_list_hr_bulk)
            #     c2 = Counter(nn_list_hr)
            #     diff = c1 - c2
            #     aa = list(diff.elements())
            #     if aa:
            #         nn_dict_bot[site_index] = [
            #             (str(site.specie), a[0], a[1], Shaper.calculate_bv(cr[str(site.specie)],
            #                                                          cr[a[0]], a[1])) for a in aa]
            #         # nn_dict_bot[site_index] = [(str(site.specie), a[0], a[1]) for a in aa]
            # nn_dict_bot = [aa[-1] for sublist in nn_dict_bot.values() for aa in sublist]
            # bbs[(nn_method, weight)][shift]['bot'] = sum(nn_dict_bot) / area
        return bbs, area

    @staticmethod
    def get_bvs(bbs, weight):
        cr = CovalentRadius().radius
        bvs_by_shift = {}
        for shift, bonds in bbs.items():
            shift_tot = 0
            for site, nn_list in bonds.items():
                specie = site.split('-')[0]
                if weight == 'BVS':
                    shift_tot += sum(
                        [Shaper.calculate_bv(cr[specie], cr[nn[0].split('-')[0]], nn[1]) for nn in nn_list])
                elif weight == 'equal':
                    shift_tot += len(nn_list)
                else:
                    raise ValueError(f'weight {weight} not recognized!')
            bvs_by_shift[shift] = shift_tot
        return bvs_by_shift

    @staticmethod
    def get_neighbors(bulk_conv, cutoff, nn_method, site_index):
        if nn_method == 'BNN':
            bnn = BrunnerNN_real(cutoff=cutoff)
            try:
                nn_list = bnn.get_nn_info(bulk_conv, site_index)
            except ValueError:
                bnn = BrunnerNN_real(cutoff=2 * cutoff)
                nn_list = bnn.get_nn_info(bulk_conv, site_index)
            nn_list = [nn['site'] for nn in nn_list]
        else:
            site = bulk_conv[site_index]
            nn_list = bulk_conv.get_neighbors(site=site, r=cutoff)
        return nn_list

    @staticmethod
    def fix_regions(struct: Union[Slab, Structure], tol: float = 0.1,
                    fix_type: str = 'z_pos') -> Union[Slab, Structure]:
        """
        Method to add site properties to the structure to fix certain ions
        from moving during a relaxation run. Mostly used to reduce computation
        time for larger structures by fixing inner layers to simulate bulk.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Main object in pymatgen to store structures.
        tol : float, optional
            Tolerance used to cluster sites into layers. The default is 0.1.
        fix_type : str, optional
            Type of region fixing to be used. For PES calculations, 'z_pos'
            is usually employed while to simulate bulk in certain slabs,
            we can fix part of the slab completely. The default is 'z_pos'.

        Raises
        ------
        Exception
            fix_type should be one of 'z_pos', 'top_half', 'bottom_half',
            'top_third', 'bottom_third', 'custom'

        Returns
        -------
        struct : pymatgen.core.structure.Structure
            Copy of the input structure with added site properties that
            adds selective dynamics properties to sites.

        """
        allowed_fix_types = ['z_pos', 'top_half', 'bottom_half', 'top_third',
                             'bottom_third']
        if fix_type not in allowed_fix_types:
            raise ValueError('Your fix_type is not in allowed_fix_types')
        layers = Shaper.get_layers(struct, tol)
        sorted_layers = sorted(layers.keys())
        num_layers = len(layers)
        fix_arr = np.asarray([[True, True, True] for _ in range(len(struct))])
        if fix_type == 'z_pos':
            fix_arr[:, 2] = False
        elif fix_type in ('top_half', 'bottom_half', 'top_third', 'bottom_third'):
            num_layers_fix = int(num_layers / 2) if fix_type.endswith('half') \
                else int(num_layers / 3)
            fixed_layers = sorted_layers[:num_layers_fix] if fix_type.startswith('bottom') \
                else sorted_layers[num_layers_fix:]
            sites = [item for sublist in [v for k, v in layers.items() if k in fixed_layers]
                     for item in sublist]
            for site in sites:
                fix_arr[site] = [False, False, False]
        struct_copy = struct.copy()
        fix_arr = fix_arr.tolist()
        struct_copy.add_site_property('selective_dynamics', fix_arr)
        return struct_copy

    @staticmethod
    def identify_slab(slab: Slab) -> dict:
        """
        Identifies the symmetry and the stoichiometry of the given slab.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Pymatgen Slab object.

        Returns
        -------
        dict
            Dictionary containing the symmetry and stoichiometry info.

        """
        sym = slab.is_symmetric()
        bulk = slab.oriented_unit_cell
        slab_formula = slab.composition.reduced_formula
        bulk_formula = bulk.composition.reduced_formula
        sto = (slab_formula == bulk_formula)
        return {'symmetric': sym, 'stoichiometric': sto}

    @staticmethod
    def generate_slabs(bulk_conv: Structure, sg_params: dict, to_file=False) -> tuple[dict, dict]:
        """
        Generates slabs with the given parameters.

        Parameters
        ----------
        bulk_conv : pymatgen.core.structure.Structure
            Conventional standard bulk structure from which to generate
            slabs from.
        sg_params : dict
            Parameters to be used in the SlabGenerator.
            Required keys:
                miller,
                slab_thick,
                vac_thick,
                max_normal_search,
                tol
            Optional keys:
                lll_reduce,
                center_slab,
                in_unit_planes,
                prim,

            Keys that concern various methods in Shaper
                resize: value -> bool, optional
                Context: Pymatgen's slab generation algorithm works by replicating the
                         oriented unit cell(OUC) a number of times to each the minimum slab size.
                         However, it does not determine the number of layers in the OUC correctly.
                         This leads to much larger slabs than one asks for, but it ensures that
                         the terminations of the top and bottom are always complementary, which
                         can be useful at times.

                         This tag determines whether to remove layers from the bottom until
                         the desired thickness is reached. It also modifies the vacuum by
                         modifying the c parameter so that we have the desired vacuum.
                         Defaults to 'False'.

                preserve_terminations: value -> bool, optional
                Context: Pymatgen sometimes generates more than one unique slab for a given miller
                         index, meaning there are distinct terminations for each one. This tag
                         determines whether to preserve the terminations when resizing slabs.
                         This is done accomplished by removing layers in 'chunks' where each chunk
                         corresponds to a portion of the oriented unit cell with height equal to
                         the distance between miller planes for the given conventional bulk structure.
                         Defaults to 'None'.

                min_thick_A: value -> float, optional
                Context: For certain structures, the layers are very close to each other (depending
                         on the tol parameter of course), where one could have a slab with many layers
                         but still thin in the physical space. This could lead to convergence issues
                         in surface energies. This tag allows the user to set a minimum slab thickness
                         in units of Angstroms, and the resizing algorithm will not remove any more layers
                         that will reduce the slab thickness below min_thick_A.
                         Defaults to 'None'.

            For more info about the description of the parameters concerning pymatgen's slabgenerator,
            please refer to the documentation of pymatgen.core.surface.SlabGenerator.
        to_file : bool, optional
            Whether to export the generated structures as VASP formatted files.
            Filenames are in the format {formula}_{miller_index}_{termination_index}.
            The default is 'False'.
        Returns
        -------
        slabs : list
            List of all the slabs consistent with the given parameters
        sg_dict: dict
            Dict with keys as miller indices and values corresponding
            pymatgen.core.surface.SlabGenerator objects.

        """
        # First, we check sg_params if it has the required keys and fill the missing
        # ones with the default values.
        # TODO: Find a better way, this is ugly.
        sg_params = check_input({'sg_params': sg_params}, ['sg_params'])['sg_params']

        max_index = sg_params.get('max_index')
        miller = sg_params.get('miller')
        if miller and not max_index:
            print(f'Generating slabs for the following miller indices: {miller}')
            miller = sg_params.get('miller')
            if isinstance(miller[0], int):
                miller = [(*miller,)]
            else:
                miller = [(*m,) for m in miller]
        elif max_index and not miller:
            miller = get_symmetrically_distinct_miller_indices(bulk_conv, max_index)
        elif max_index and miller:
            raise ValueError('You cannot specify both max_index and miller parameters.')
        else:
            raise ValueError('You must specify either max_index or miller parameters.')

        tol = sg_params.get('tol')
        resize = sg_params.get('resize')
        symmetrize = sg_params.get('symmetrize')
        match_ouc_lattice = sg_params.get('match_ouc_lattice')
        calculate_bonds = sg_params.get('calculate_bonds')
        sg_dict = {}
        slabs_dict = {}
        for m in miller:
            # we need some parameters to use in SlabGenerator, so we extract those
            # from the input sg_params and put them in pmg_sg_params
            pmg_sg_params = get_pmg_sg_params(bulk_conv=bulk_conv, miller=m, sg_params=sg_params)

            # we first try a SlabGenerator with the given sg_params, if things go as
            # expected we proceed with this
            sg = SlabGenerator(**pmg_sg_params)

            d_hkl, pmg_layer_size = Shaper.get_pmg_layer_size(bulk_conv=bulk_conv,
                                                              miller=m,
                                                              sg=sg,
                                                              tol=tol)

            min_thick_A = sg_params['min_thick_A']
            # if there is a min_thick_A key in sg_params, we have to ensure that the slabs we initially
            # generate (before resizing) have thicknesses greater than this value. For this, we modify
            # the corresponding min_slab_size parameter in pmg_sg_params by calculating the number of layers
            # needed to reach min_thick_A.
            if min_thick_A:
                final_layer_spacing = Shaper.get_layer_spacings(sg.oriented_unit_cell, tol)[-1]
                min_slab_size = max(np.ceil((min_thick_A + final_layer_spacing) / d_hkl), sg_params['slab_thick'])
                pmg_sg_params['min_slab_size'] = min_slab_size
                sg = SlabGenerator(**pmg_sg_params)

            slabs = sg.get_slabs(ftol=tol, symmetrize=symmetrize)
            if not slabs:
                continue

            if match_ouc_lattice:
                # we check if we can get an oriented unit cell with the same lateral lattice
                # parameters and gamma angle as the slab, for consistency with brillouin zone
                # sampling
                ouc = Shaper.get_matching_ouc(slabs[0])
                if ouc:
                    param_modified = False
                else:
                    # if no such ouc exists, we turn off primitive and lll_reduce, since
                    # only when either one or both of these are True, we have issues with
                    # finding matching ouc
                    print('OUC matching failed with the given sg_params, modifying primitive and lll_reduce..')
                    pmg_sg_params.update({'primitive': False,
                                          'lll_reduce': False})
                    sg_modified = SlabGenerator(**pmg_sg_params)
                    slabs_modified = sg_modified.get_slabs(ftol=tol, symmetrize=symmetrize)
                    # we set a flag to show that sg params are modified, and we will assign
                    # this as an attribute to the Slab objects, so that we can tell if the slabs
                    # generated result from a modified sg
                    ouc = Shaper.get_matching_ouc(slabs_modified[0])
                    if ouc:
                        param_modified = True
                        slabs = slabs_modified
                    else:
                        # if we still don't have a matching ouc, which should not happen
                        # we print and a non-matching ouc is used instead.
                        print(f'Matching oriented unit cell cannot be found. Your reference energies'
                              f'might not be suitable for a surface energy convergence scheme.')
                        param_modified = False

                # we change the oriented unit cells of slabs to the matching ouc that we find.
                if ouc:
                    for slab in slabs:
                        slab.oriented_unit_cell = ouc

                # we assign attributes to slabs if they result from a modified sg, and this is done after resizing
                # because as soon as a pymatgen structure is copied (such is the case in Shaper.resize()), it loses
                # all attributes not defined in the copy method. Since param_modified is such an attribute, we need
                # to add it after resizing the slabs.
            else:
                param_modified = False

            for slab in slabs:
                slab.param_modified = param_modified

            # resize flag is used generate slabs with user-defined thicknesses in number of layers.
            # This is done by removing layers from the bottom of the pymatgen generated slabs since
            # they are usually thicker than one would expect.
            slab_thick, vac_thick = sg_params['slab_thick'], sg_params['vac_thick']
            if resize:
                # if we want to preserve terminations, we must remove layers in chunks
                # and these chunk sizes are determined by pmg_layer_size
                preserve_terminations = sg_params.get('preserve_terminations')
                chunk_size = pmg_layer_size if preserve_terminations else 1

                slabs = [Shaper.resize(slab, slab_thick, vac_thick, tol=tol,
                                       chunk_size=chunk_size, min_thick_A=min_thick_A)
                         for slab in slabs]
            else:
                # TODO: Remove this once resize is confirmed working. Only here to generate
                # TODO: slabs with sizes comparable to the initial slab_thick with pymatgen
                # slab_thick_pmg = np.ceil(slab_thick / pmg_layer_size)
                # pmg_sg_params.update({'min_slab_size': slab_thick_pmg})
                # sg = SlabGenerator(**pmg_sg_params)
                # slabs = sg.get_slabs(ftol=tol, symmetrize=symmetrize)
                # THIS NOW JUST RESIZES THE VACUUM AND LEAVES THE SLABS UNTOUCHED
                slabs = [Shaper.resize(slab, vacuum_thickness=vac_thick, tol=tol) for slab in slabs]

            # check if slabs list is empty, which only happens when symmetrize is True
            # and the sg can not find symmetric slabs.
            try:
                slab = slabs[0]
            except IndexError:
                print(f'Symmetric slabs could not be generated for {m} orientation. Increasing slab_thick'
                      ' may or may not solve this issue.')
                continue

            for slab in slabs:
                slab.pmg_layer_size = pmg_layer_size

            # TODO: Remove this once we decide where to filter out large slabs
            # max_nsites = sg_params.get('max_nsites', None)
            # if max_nsites:
            #     slabs = [slab for slab in slabs if slab.num_sites <= max_nsites]
            # try:
            #     slab = slabs[0]
            # except IndexError:
            #     print(f'No slabs could be generated for {m} orientation because all'
            #           f' slabs have more than {max_nsites} sites.')
            #     continue

            # we assign energies (bond valence sums of broken bonds) to each slab, in this case
            # unique terminations, and also the pymatgen layer size, which is the number of layers
            # we can safely remove at once while preserving terminations.

            if calculate_bonds:
                nn_method = sg_params.get('nn_method')
                weight = sg_params.get('weight')
                max_bl = max([a[0] for a in list(Shaper.get_bonds(bulk_conv).values())])
                bbs, area = Shaper.bonds_by_shift_new(sg=sg,
                                                      bulk_conv=bulk_conv,
                                                      nn_method=nn_method,
                                                      tol=tol,
                                                      cutoff=max_bl,
                                                      edge_tol=99)
                bvs = Shaper.get_bvs(bbs, weight=weight)
                for slab in slabs:
                    slab.energy = {'broken_bonds': bbs[np.round(slab.shift, 4)],
                                   'bvs_per_area': bvs[np.round(slab.shift, 4)] / area,
                                   'area': area}

            if to_file:
                formula = slabs[0].composition.reduced_formula
                for index, slab in enumerate(slabs):
                    hkl = ''.join([str(i) for i in slab.miller_index])
                    area = np.round(slab.surface_area, 2)
                    slab.to(f'{formula}_{hkl}_{area}_{index}.vasp', 'poscar')

            slabs_dict[m] = slabs
            sg_dict[m] = sg

        return slabs_dict, sg_dict

    @staticmethod
    def get_pmg_layer_size(bulk_conv, miller, sg, tol):
        # since the terminations repeat every d_hkl distance in c direction,
        # the distance between miller planes, we need to figure out how many
        # layers this d_hkl portion corresponds to in order to preserve terminations
        d_hkl = bulk_conv.lattice.d_hkl(miller)
        try:
            ouc_layers = len(Shaper.get_layers(sg.oriented_unit_cell, tol))
        except ValueError:
            ouc_layers = 1
        ouc_height = Shaper.get_proj_height(sg.oriented_unit_cell)
        # we calculate how many layers pymatgen considers a single layer here
        pmg_layer_size = int(ouc_layers / round(ouc_height / d_hkl))
        return d_hkl, pmg_layer_size

    @staticmethod
    def get_constrained_ouc(slab: Slab) -> Structure:
        """
        Finds the constrained oriented unit cell of a Slab object. The constraints
        are the a and b parameters of the slab, along with the gamma angle.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Pymatgen Slab object whose oriented unit cell we want to constrain.

        Returns
        -------
        ouc : pymatgen.core.structure.Structure
            Constrained oriented unit cell of the input Slab.

        """
        constraints = {'a': slab.lattice.a,
                       'b': slab.lattice.b,
                       'gamma': slab.lattice.gamma}
        ouc = slab.oriented_unit_cell.get_primitive_structure(constrain_latt=constraints)
        return ouc

    @staticmethod
    def get_matching_ouc(slab: Slab) -> Union[Structure, None]:
        """
        Given a slab, finds an oriented unit cell that matches the lateral lattice parameters
        and the gamma angle of the slab. Useful for constructing an oriented unit cell to be
        used as a reference structure for surface energy calculations.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Pymatgen Slab object whose oriented unit cell we want to constrain.

        Returns
        -------
        ouc : pymatgen.core.structure.Structure
            Constrained oriented unit cell of the input Slab.

        """
        # applying LLL reduction on a structure sometimes changes the orders of the lattice
        # parameters and hence, the ordering of the lattice vectors. In order to have a
        # consistent sampling of Brillouin zone between the slab and the oriented unit cell
        # we rotate the oriented unit cell to have the same orientation as the slab.
        trans = {0: ((0, 0, 1), (0, 1, 0), (1, 0, 0)),
                 1: ((1, 0, 0), (0, 0, 1), (0, 1, 0))}
        ouc = slab.oriented_unit_cell
        # we first check if the preset OUC matches
        lattice_match = Shaper._check_lattice_match(slab.lattice, ouc.lattice)
        # angle_check = ouc.lattice.angles[3 - sum(indices)] == slab.lattice.gamma
        if not lattice_match:
            ouc = ouc.copy(sanitize=True)
            lattice_match = Shaper._check_lattice_match(slab.lattice, ouc.lattice)
            if not lattice_match:
                ouc = Shaper.get_constrained_ouc(slab)
                lattice_match = Shaper._check_lattice_match(slab.lattice, ouc.lattice)
                if not lattice_match:
                    return None

        # since the whole reason we do this is to match the kpoint sampling of the slab
        # and the oriented unit cell in the plane parallel to the surface, we need to align
        # the matched lattice vectors, and that's why we perform this SupercellTransformation
        # to make sure we transform the OUC to have the same base vectors as the slab
        if lattice_match != 2:
            st = SupercellTransformation(trans[lattice_match])
            ouc = st.apply_transformation(ouc)
        return ouc

    @staticmethod
    def _check_lattice_match(lattice1: Lattice, lattice2: Lattice) -> Union[int, None]:
        """
        checks if lattice1 has the same base as lattice2, ignoring orientations.
        lattice1 is considered to be the slab in this case. returns the index of
        the third vector that is not shared between the lattices.

            Parameters
            ----------
            lattice1 : pymatgen.core.lattice.Lattice
                Pymatgen Lattice object that we want to use as a reference.
            lattice2 : pymatgen.core.lattice.Lattice
                Pymatgen Lattice object that we want to compare against the reference.

            Returns
            -------
            angle_index (int) or None
            if there is a match, we return the index of the angle between the matched
            base vectors. This is enough information to go further. If no match is found,
            we return None

            """
        matches_ab = get_subset_indices(lattice1.abc[:2], lattice2.abc)
        if not matches_ab:
            return None
        else:
            for match_ab in matches_ab:
                angle_index = 3 - sum(match_ab)
                check = np.isclose(lattice1.gamma, lattice2.angles[angle_index])
                if check:
                    return angle_index
            return None
