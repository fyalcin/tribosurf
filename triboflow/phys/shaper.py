#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 01:05:05 2021

Class and methods that deal with the general Shaper of structures.

The module contains:

    ** Shaper **:
        General class to examine layers, bonds, lattice parameters,
        and reconstructs slabs with desired transformations.
        It includes the following methods:
            - _get_layer_spacings
            - reconstruct_slab
            - get_surface_normal
            - _get_hkl_projection
            - _get_layers
            - _remove_layers

    Functions:
    - attr_to_dict

    Author: Fırat Yalçın

"""

__author__ = 'Fırat Yalçın'
__contact__ = 'firat.yalcin@univie.ac.at'
__date__ = 'April 21st, 2021'

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import center_slab, Slab, SlabGenerator, get_symmetrically_distinct_miller_indices
from pymatgen.analysis.local_env import BrunnerNN_real
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.io.cif import CifParser
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from collections import defaultdict
import numpy as np
import itertools


def range_diff(r1, r2):
    s1, e1 = r1
    s2, e2 = r2
    endpoints = sorted((s1, s2, e1, e2))
    result = []
    if endpoints[0] == s1:
        result.append((endpoints[0], endpoints[1]))
    if endpoints[3] == e1:
        result.append((endpoints[2], endpoints[3]))
    return result


def multirange_diff(r1_list, r2_list):
    for r2 in r2_list:
        r1_list = list(itertools.chain(*[range_diff(r1, r2) for r1 in r1_list]))
    return r1_list


def attr_to_dict(obj, attrs):
    attr_dict = {attr: getattr(obj, attr, None) for attr in attrs}
    return attr_dict


class Shaper():

    @staticmethod
    def _get_layer_spacings(struct, tol=0.1):
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
    
        Returns
        -------
        list
            list of floats representing the projected distances between layers
            along the surface normal direction in angstroms
    
        """

        # Layer info that contains the c-coordinates and sites
        layers = Shaper._get_layers(struct, tol)

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
        proj_height = Shaper._get_proj_height(struct)

        return np.round([spacing * proj_height for spacing in d], 10)

    @staticmethod
    def _get_proj_height(struct, region='cell', min_vac=4.0):
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

        proj_height = Shaper._get_hkl_projection(struct.lattice.matrix[2], struct)
        if region == "cell":
            return proj_height
        elif region == "slab" or region == "vacuum":
            spacings = Shaper._get_layer_spacings(struct)
            slab_height = sum([s for s in spacings if s < min_vac])
            return slab_height if region == "slab" else proj_height - slab_height
        else:
            raise ValueError('Region must be one of "cell", "vacuum", or "slab"')

    @staticmethod
    def reconstruct(struct, struct_thickness, vacuum_thickness, tol=0.1, minimize_bv=True,
                    center=True, **kwargs):
        """
        Reconstruct the input slab with the desired slab thickness in
        number of layers and the vacuum region in Angstroms. All the attributes
        of sites are preserved by the reconstruction.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Structure object that is to be reconstructed. Input object
            is not modified with this method.
        slab_thickness : int
            Desired slab thickness in number of layers. Layers will
            be removed from the bottom until the desired thickness is
            reached.
        vacuum_thickness : float
            Desired vacuum region thickness in Angstroms. Lattice
            parameters are modified in order to get the correct vacuum.
        minimize_bv : bool, optional
            Whether to minimize the bond valence sum of broken bonds when
            removing layers. The default is True.
        center : bool
            Whether to center the reconstructed slab between the vacuum region.
            The default is True.

        Returns
        -------
        reconstructed_struct : pymatgen.core.structure.Structure
            Reconstructed structure with the desired parameters.

        """
        # Input slab is first centered for the cases where the slab spills
        # outside the box from the top and the bottom
        struct_centered = center_slab(struct.copy(sanitize=True))
        initial_thickness = Shaper._get_proj_height(struct_centered, 'slab')
        # Layers (containing sites) are removed from the bottom until
        # the desired slab_thickness is reached
        if minimize_bv:
            bbs = kwargs['bbs']
            bvs, indices = np.unique(list(bbs.values()), return_index=True)
            periodicity = len(bvs)
            spacings = [spacing for spacing in Shaper._get_layer_spacings(struct_centered, tol) if spacing < 4.0]
            num_layers = len(Shaper._get_layers(struct_centered, tol))
            layers_to_remove = int(periodicity * np.floor((num_layers
                                                           - struct_thickness) / periodicity))
            while initial_thickness - sum(spacings[:layers_to_remove]) < 10:
                if layers_to_remove < 0:
                    raise ValueError('you must choose a bigger slab')
                layers_to_remove -= periodicity
            struct_resized = Shaper._remove_layers(struct_centered, layers_to_remove,
                                                   tol=tol, method='layers')
            if struct_resized.composition.reduced_composition != struct_centered.composition.reduced_composition:
                print('hey')
        else:
            struct_resized = Shaper._remove_layers(struct_centered, struct_thickness, tol=tol)

        # Vacuum region is modified to the desired thickness
        reconstructed_struct = Shaper._modify_vacuum(struct_resized, vacuum_thickness)
        if reconstructed_struct.composition.reduced_composition != struct_centered.composition.reduced_composition:
            print('hey')

        return reconstructed_struct

    @staticmethod
    def _modify_vacuum(struct, vac_thick, method='to_value', center=True):
        """
        Method to modify the vacuum region in a structure.
    
        Parameters
        ----------
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
            # Necessary slab attributes to reconstruct the Slab
            attrs = ["species", "miller_index", "oriented_unit_cell",
                     "shift", "scale_factor", "reorient_lattice",
                     "reconstruction", "site_properties", "energy"]
            struct_params = attr_to_dict(struct, attrs)
            out_object = Slab
        else:
            # Necessary structure attributes to reconstruct the Structure
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
        proj_height = Shaper._get_proj_height(struct)

        # 'c' parameter of the Lattice is modified to adjust vacuum
        # to the desired thickness
        if method == 'to_value':
            initial_vac = Shaper._get_proj_height(struct, 'vacuum')
            lat_params['c'] += (vac_thick - initial_vac) * lat_params['c'] / proj_height
        elif method == 'by_value':
            lat_params['c'] += vac_thick * lat_params['c'] / proj_height

        new_lat = Lattice.from_parameters(**lat_params)

        modified_struct = center_slab(out_object(new_lat, **struct_params)) \
            if center else out_object(new_lat, **struct_params)

        return modified_struct

    @staticmethod
    def _get_hkl_projection(vector, struct):
        """
        Simple method to calculate the norm of theprojection of a vector
        along the hkl direction which is assumed to be normal to the plane
        formed by the first two lattice vectors of the passed structure.
        Useful for structures where the third lattice vector is not in
        the same direction as the surface normal.
    
        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Main object in pymatgen to store structures.
    
        Returns
        -------
        float
            Projected height of the given structure in the direction
            that is normal to the x-y plane
    
        """
        latvec = struct.lattice.matrix
        normal = np.cross(latvec[0], latvec[1])
        normal /= np.linalg.norm(normal)
        return np.dot(vector, normal)

    @staticmethod
    def _get_layers(struct, tol=0.1):
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
                cdist = frac_coords[i][2] - frac_coords[j][2]
                # cdist = abs(cdist - round(cdist)) * proj_height
                cdist = abs(cdist - round(cdist)) * struct.lattice.abc[2]
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

        # for each layer, find sites that belong to it and assign an
        # average c-value for the layer
        layers = {sum([struct.frac_coords[i][2] - np.floor(struct.frac_coords[i][2])
                       for i in v]) / len(v): v for k, v in layers.items()}
        return layers

    @staticmethod
    def _remove_layers(slab, num_layers, tol=0.1, method='target', position='bottom',
                       center=True):
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
        layers = Shaper._get_layers(slab, tol)
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
        indices_list = [layers[c_coord] for c_coord in to_remove]
        flat_list = [item for sublist in indices_list for item in sublist]
        slab_copy = slab.copy()
        slab_copy.remove_sites(flat_list)
        return center_slab(slab_copy) if center else slab_copy

    @staticmethod
    def _get_average_layer_spacing(slab, tol=0.1, vacuum_treshold=6.0):
        """
        Compute the average distance between the slabs layers disregarding the
        vacuum region.
    
        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Standard pymatgen Slab object.
        vacuum_treshold : float, optional
            Regions larger than this will be treated as vacuum and will not be
            treated as an interlayer spacing. The default is 6.0
        Returns
        -------
        av_spacing : float
            Average layer spacing
    
        """
        spacings = Shaper._get_layer_spacings(slab, tol)
        spacings_no_vac = np.delete(spacings,
                                    np.where(spacings >= vacuum_treshold))
        av_spacing = np.mean(spacings_no_vac)
        return av_spacing

    @staticmethod
    def _get_bonds(struct, method='covalent_radii', dtol=0.20, wtol=0.15):
        """
        Finds all unique bonds in the structure and orders them by bond strength
        using bond valance method and with the assumption that the ideal bond length
        = CovalentRadius(site1) + CovalentRadius(site2)

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Conventional standart structure that is used to generate the slabs.

        method : string, optional
            Method used to calculate the bond valence parameters
            - 'covalent_radii' : sets the 'ideal' bond length R_0 as the sum of
                the covalent radii of bonded atoms
            - 'BVparams' : uses a .cif file of a list of fitted bond valence
                parameters from various sources.

        dtol : float, optional
            Added tolerance to form a bond for the dictionary passed to the
            slab generation algorithm.

        wtol : float, optional
            Added tolerance to eliminate bonds by their weights calculated by
            exp((R_0 - R_i/b)) where R_0 is the 'ideal' bond length, R_i is
            the observed bond length in the structure, and b is an empirical
            constant roughly 0.37 Angstroms.

        Returns
        -------
        dict : Collection of bonds that has a 'weight' within a delta of the highest
        weight.
        """
        # struct = struct.get_primitive_structure()
        BNN = BrunnerNN_real(cutoff=max(struct.lattice.abc))
        species, indices = np.unique([str(x) for x in struct.species],
                                     return_index=True)
        bonds = {}
        wmax = 0
        for i, site_index in enumerate(indices):
            sp1 = species[i]
            for neighbor in BNN.get_nn_info(struct, site_index):
                sp2 = str(struct.species[neighbor['site_index']])
                dist = np.linalg.norm(struct[site_index].coords
                                      - neighbor['site'].coords)
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
                    bonds[(sp1, sp2)] = (dist + dtol, w)
        # bonds = {k: v[0] for k, v in bonds.items() if abs(v[1]-wmax)/wmax <= wtol}
        return bonds

    @staticmethod
    def _get_c_ranges(struct, nn_method='all'):
        """
        Calculates all the bonds in the given structure

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Pymatgen Structure object.
        nn_method: str, optional
            Nearest-neighbor algorithm to be used. Currently supports
            'all' and 'BNN'. For more info, check out BrunnerNN documentation.

        Returns
        -------
        c_ranges : list
            List with elements describing every bond in the structure, with
            the endpoints of the bond and the bond valence values for each bond.

        """
        cr = CovalentRadius().radius
        if nn_method == 'all':
            nn_list = struct.get_all_neighbors(max(struct.lattice.abc))
        elif nn_method == 'BNN':
            BNN = BrunnerNN_real(cutoff=max(struct.lattice.abc))
            nn_list = BNN.get_all_nn_info(struct)
        c_ranges = []
        for s_index, site in enumerate(struct):
            for nn in nn_list[s_index]:
                c_range = np.round(sorted([site.frac_coords[2], nn.frac_coords[2]]), 3)
                if c_range[0] != c_range[1]:
                    nn_site_index = nn.index
                    sp1 = str(struct.species[s_index])
                    sp2 = str(struct.species[nn_site_index])
                    dist = nn.nn_distance
                    bv = Shaper._get_bv(cr[sp1], cr[sp2], dist)
                    bv = ((sp1, s_index), (sp2, nn_site_index), dist, bv)
                    if c_range[0] < 0:
                        c_ranges.append((0, c_range[1], bv))
                        c_ranges.append((c_range[0] + 1, 1, bv))
                    elif c_range[1] > 1:
                        c_ranges.append((c_range[0], 1, bv))
                        c_ranges.append((0, c_range[1] - 1, bv))
                    else:
                        c_ranges.append((c_range[0], c_range[1], bv))
        return c_ranges

    @staticmethod
    def _get_bv(r1, r2, bond_dist):
        b = 0.37
        R_0 = r1 + r2
        return np.exp((R_0 - bond_dist) / b)
    
    @staticmethod
    def get_surface_area(struct):
        mat = struct.lattice.matrix
        return np.linalg.norm(np.cross(mat[0], mat[1]))

    @staticmethod
    def _bonds_by_shift(SG, nn_method='all', tol=0.1):
        """
        Calculates the bond valence sums of the broken bonds corresponding to
        all the possible shifts

        Parameters
        ----------
        SG : pymatgen.core.surface.SlabGenerator
            Pymatgen SlabGenerator object to extract the possible shifts in a
            specific orientation.
        tol : float, optional
            Tolerance value used in the layering of sites in units of Angstroms.
            The default is 0.1.

        Returns
        -------
        bbs : dict
            Dictionary with keys as shifts and values as the bond valence sum of the
            broken bonds at each shift, scaled by the area of the x-y plane.

        """
        ouc = SG.oriented_unit_cell
        area = Shaper.get_surface_area(ouc)
        c_ranges = Shaper._get_c_ranges(ouc, nn_method)
        shifts = np.round(SG._calculate_possible_shifts(tol=tol), 4)
        bbs = {}
        for shift in shifts:
            bbs[shift] = 0
            for c_range in c_ranges:
                if c_range[0] < shift < c_range[1]:
                    bbs[shift] += c_range[2][3]
            bbs[shift] = np.round(bbs[shift], 4)/area
        return bbs

    @staticmethod
    def fix_regions(struct, tol=0.1, fix_type='z_pos'):
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
            is usually comployed while to simulate bulk in certain slabs,
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
        layers = Shaper._get_layers(struct, tol)
        sorted_layers = sorted(layers.keys())
        num_layers = len(layers)
        fix_arr = np.asarray([[True, True, True] for i in range(len(struct))])
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
        struct_copy.add_site_property('selective_dynamics', fix_arr)
        return struct_copy

    @staticmethod
    def identify_slab(slab):
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
        bulk_ref = slab.oriented_unit_cell
        slab_comp = slab.composition.reduced_composition
        bulk_ref_comp = bulk_ref.composition.reduced_composition
        sto = (slab_comp == bulk_ref_comp)
        return {'symmetric': sym, 'stoichiometric': sto}

    @staticmethod
    def generate_slabs(bulk_conv, sg_params, reconstruct=True, to_file=False):
        """
        Generates slabs with the given parameters.

        Parameters
        ----------
        bulk_conv : pymatgen.core.structure.Structure
            Conventional standard bulk structure from which to generate
            slabs from.
        sg_params : dict
            Parameters to be used in the SlabGenerator.
            Required keys are:
                miller,
                slab_thick,
                vac_thick,
                max_normal_search,
                tol
            Optional keys are:
                lll_reduce,
                center_slab,
                in_unit_planes,
                prim,
                minimize_bv

            For more info about the description of these parameters,
            refer to the documentation of pymatgen.core.surface.SlabGenerator.
        reconstruct : bool, optional
            Context: Pymatgen's slab generation algorithm works by replicating the
            oriented unit cell(OUC) a number of times to each the minimum slab size.
            However, it does not determine the number of layers in the OUC correctly.
            This leads to much larger slabs than one asks for, but it ensures that
            the terminations of the top and bottom are always complementary, which
            can be useful at times.

            This tag determines whether to remove layers from the bottom until
            the desired thickness is reached. It also modifies the vacuum by
            modifying the c parameter so that we have the desired vacuum.
            The default is 'True'.

        to_file : bool, optional
            Whether to export the generated structures as VASP formatted files.
            Filenames are in thr format {formula}_{miller_index}_{termination_index}.
            The default is 'False'.
        Returns
        -------
        slabs : list
            List of all the slabs (unique terminations) generated with the given
            parameters
        SG : pymatgen.core.surface.SlabGenerator
            Pymatgen SlabGenerator object used to generate the slabs.

        """
        max_index = sg_params.get('max_index')
        if max_index:
            miller = get_symmetrically_distinct_miller_indices(bulk_conv, max_index)
        else:
            miller = sg_params.get('miller')
            if isinstance(miller[0], int):
                miller = [(*miller,)]
            else:
                miller = [(*m,) for m in miller]
        slab_thick = sg_params.get('slab_thick')
        vac_thick = sg_params.get('vac_thick')
        minimize_bv = sg_params.get('minimize_bv')
        mns = sg_params.get('max_normal_search')

        SG_dict = {}
        slabs_list = []
        for m in miller:
            max_normal_search = max([abs(i) for i in m]) if mns == 'max' else mns
            SG = SlabGenerator(initial_structure=bulk_conv,
                               miller_index=m,
                               min_slab_size=slab_thick,
                               min_vacuum_size=vac_thick,
                               lll_reduce=sg_params.get('lll_reduce', True),
                               center_slab=sg_params.get('center_slab', True),
                               in_unit_planes=sg_params.get('in_unit_planes', True),
                               primitive=sg_params.get('prim', True),
                               max_normal_search=max_normal_search,
                               reorient_lattice=True)

            nn_method = 'all'
            tol = sg_params.get('tol')
            bbs = Shaper._bonds_by_shift(SG, nn_method, tol)
            slabs = SG.get_slabs(tol=tol, symmetrize=sg_params.get('symmetrize', False))

            if reconstruct:
                slabs = [Shaper.reconstruct(slab, slab_thick, vac_thick, tol=tol,
                                            minimize_bv=minimize_bv, bbs=bbs)
                         for slab in slabs]

            for slab in slabs:
                slab.energy = bbs[np.round(slab.shift, 4)]

            if to_file:
                formula = slabs[0].composition.reduced_formula
                for index, slab in enumerate(slabs):
                    hkl = ''.join([str(i) for i in slab.miller_index])
                    area = np.round(slab.surface_area, 2)

                    slab.to('poscar', f'{formula}_{hkl}_{area}_{index}.vasp')

            slabs_list += slabs
            SG_dict[m] = SG

        return slabs_list, SG_dict


    @staticmethod
    def get_constrained_ouc(slab):
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
