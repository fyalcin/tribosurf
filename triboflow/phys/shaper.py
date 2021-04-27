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
from pymatgen.core.surface import center_slab, Slab, get_slab_regions
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from collections import defaultdict
import numpy as np
import itertools


def attr_to_dict(obj, attrs):
    attr_dict = {attr: getattr(obj, attr, None) for attr in attrs}
    return attr_dict


class Shaper():

    @staticmethod
    def _get_layer_spacings(struct, tol=0.1):
        """
        Simple method to calculate the thickness of the vacuum region.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Pymatgen object to store slabs. Note that the slab
            should be oriented in such a way that the surface should
            be parallel to the plane that the first 2 lattice vectors lie on.

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
        d = [x - layers_c[i-1] for i, x in enumerate(layers_c)]

        # When a periodic boundary is passed, layers wrap over and we get a
        # negative spacing, to correct, we add 1 to negative spacing values
        d = [s + int(s < 0) for s in d]

        # For slabs with the third lattice vector not along miller
        # direction, we need the projected height to also project the vacuum
        # height
        proj_height = Shaper._get_hkl_projection(struct.lattice.matrix[2], struct)

        return np.round([spacing*proj_height for spacing in d], 10)

    @staticmethod
    def _get_proj_height(struct, region='cell'):
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
            slab_regions = get_slab_regions(struct, min(struct.lattice.abc))
            slab_height = sum([proj_height*(reg[1]-reg[0]) for reg in slab_regions])
            return slab_height if region == "slab" else proj_height - slab_height
        else:
            raise ValueError('Region must be one of "cell", "vacuum", or "slab"')
        return proj_height

    @staticmethod
    def _identify_regions(struct, min_vac=4.0):
        """
        Internal method to identify regions in a given structure.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Main object in pymatgen to store structures.
        min_vac : float, optional
            Minimum thickness to identify a region as 'vacuum'.
            The default is 4.0.

        Returns
        -------
        regions : dict
            Simple dictionary with keys as regions 'slab' and 'vacuum' and values
            as the respective projected thickness values in angstroms.

        """
        frac_coords_c = struct.frac_coords[:, 2]
        struct_cp = struct.copy()
        # To avoid issues with identifying vacuum region at the periodic boundary,
        # we translate the structure in -c direction by the lowest c coordinate of
        # its constituent sites. This way, we always start with a slab region.
        struct_cp.translate_sites(list(range(len(struct_cp))), [0, 0, -min(frac_coords_c)])

        spacings = Shaper._get_layer_spacings(struct_cp)
        vac_region = sum([s for s in spacings if s > min_vac])
        slab_region = sum(spacings) - vac_region
        regions = {"vacuum": vac_region, "slab": slab_region}
        return regions

    @staticmethod
    def reconstruct(struct, struct_thickness, vacuum_thickness, center=True):
        """
        Reconstruct the input slab with the desired slab thickness in
        number of layers and the vacuum region in Angstroms. All the attributes
        of sites are preserved by the reconstruction.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Slab object that is to be reconstructed. Input object
            is not modified with this method.
        slab_thickness : int
            Desired slab thickness in number of layers. Layers will
            be removed from the bottom until the desired thickness is
            reached.
        vacuum_thickness : float
            Desired vacuum region thickness in Angstroms. Lattice
            parameters are modified in order to get the correct vacuum.
        center : bool
            Whether to center the reconstructed slab between the vacuum region.
            The default is True.

        Returns
        -------
        reconstructed_slab : pymatgen.core.surface.Slab
            Reconstructed slab with the desired parameters.

        """
        # Input slab is first centered for the cases where the slab spills
        # outside the box from the top and the bottom
        struct_centered = center_slab(struct.copy(sanitize=True))

        # Layers (containing sites) are removed from the bottom until
        # the desired slab_thickness is reached
        struct_resized = Shaper._remove_layers(struct_centered, struct_thickness)

        # Check if a Slab or Structure is passed and process accordingly
        if 'miller_index' in vars(struct_resized):
            # Necessary slab attributes to reconstruct the Slab
            attrs = ["species", "miller_index", "oriented_unit_cell",
                     "shift", "scale_factor", "reorient_lattice",
                     "reconstruction", "site_properties", "energy"]
            struct_params = attr_to_dict(struct_resized, attrs)
            out_object = Slab
        else:
            # Necessary structure attributed to reconstruct the Structure
            attrs = ["species", "site_properties"]
            struct_params = attr_to_dict(struct_resized, attrs)
            out_object = Structure

        # To avoid issues with fractional coordinates when scaling vacuum,
        # cartesian coordinates are used
        corrected_params = {'coords': struct_resized.cart_coords,
                            'coords_are_cartesian': True}
        struct_params.update(corrected_params)

        # Initial vacuum region is calculated
        initial_vacuum = Shaper._identify_regions(struct_resized).get('vacuum')

        # Lattice parameters are generated in order to be modified
        lat_attrs = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        lat_params = attr_to_dict(struct_resized.lattice, lat_attrs)

        latvec = struct_resized.lattice.matrix
        proj_height = Shaper._get_hkl_projection(latvec[2], struct_resized)

        # 'c' parameter of the Lattice is modified to adjust vacuum
        # to the desired thickness
        lat_params['c'] += (vacuum_thickness-initial_vacuum)*lat_params['c']/proj_height
        new_lat = Lattice.from_parameters(**lat_params)

        # Reconstructed slab is generated from the resized slab parameters
        # and modified Lattice, then centered again due to being off centered
        # after layer removal and vacuum modification

        reconstructed_struct = center_slab(out_object(new_lat, **struct_params)) \
            if center else out_object(new_lat, **struct_params)

        return reconstructed_struct

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
            Tolerance used to cluster sites into layers. The default is 0.1.

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
                       for i in v])/len(v): v for k, v in layers.items()}
        return layers

    @staticmethod
    def _remove_layers(slab, num_layers, method='target', position='bottom'):
        """
        Removes layers from the bottom of the slab while updating the number
        of bonds broken in the meantime.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Standard pymatgen Slab object.
        num_layers : int
            Number of layers to remove from the structure
        method : str
            Whether to remove num_layers or remove layers until the
            structure has num_layers number of layers in total.
            Options are 'target' and 'layers'. The default is 'target'.
        position : string, optional
            Side on which the sites should be removed.
            Available options are 'top' and 'bottom'. The default is 'bottom'.

        Returns
        -------
        slab_copy : pymatgen.core.surface.Slab
            Copy of the input Slab structure with layers removed.

        """
        layers = Shaper._get_layers(slab)
        if num_layers > len(layers):
            raise ValueError('Number of layers to remove/target can\'t exceed \
                             the number of layers in the given slab.')
        c_coords = sorted(layers.keys())
        if method == "layers":
            to_remove = c_coords[:num_layers] if position == "bottom" \
                        else c_coords[num_layers:]
        elif method == "target":
            to_remove = c_coords[:len(c_coords)-num_layers] if position == "bottom" \
                else c_coords[len(c_coords)-num_layers:]
        indices_list = [layers[c_coord] for c_coord in to_remove]
        flat_list = [item for sublist in indices_list for item in sublist]
        slab_copy = slab.copy()
        slab_copy.remove_sites(flat_list)
        return slab_copy
