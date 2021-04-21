#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 01:05:05 2021

Class and methods that deal with the general geometry of structures.

The module contains:

    ** Geometry **:
        General class to examine layers, bonds, lattice parameters, 
        and reconstructs slabs with desired transformations.
        It includes the following methods:
            - get_vacuum
            - reconstruct_slab
            - get_surface_normal
            - get_proj_height
            - get_layers
            - remove_layers
    
    Functions:
    - attr_to_dict

    Author: Fırat Yalçın

"""

__author__ = 'Fırat Yalçın'
__contact__ = 'firat.yalcin@univie.ac.at'
__date__ = 'April 21st, 2021'

from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import center_slab, Slab
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from collections import defaultdict
import numpy as np
import itertools
        
def attr_to_dict(obj, attrs):
    attr_dict = {attr : getattr(obj, attr, None) for attr in attrs}
    return attr_dict

class Geometry():
    
    @staticmethod
    def get_vacuum(slab):
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
        float
            Thickness of the vacuum region in Angstroms.

        """
        slab_cp = slab.copy()
        # Layer info that contains the c-coordinates and sites
        layers = Geometry.get_layers(slab_cp)
        # Only the c-coordinates of the layers are needed
        layers_c = sorted(layers.keys())
        # Spacing between consecutive layers are calculated
        d = [x - layers_c[i-1] for i,x in enumerate(layers_c)]
        # When a periodic boundary is passed, layers wrap over and
        # we get a negative spacing, to correct, we add 1 to negative 
        # spacing values
        d = [s + int(s<0) for s in d]
        # For slabs with the third lattice vector not along miller
        # direction, we need the projected height to also project the vacuum
        # height
        proj_height = Geometry.get_proj_height(slab)
        return np.round(max(d)*proj_height,4)
    
    @staticmethod
    def reconstruct_slab(slab, slab_thickness, vacuum_thickness):
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

        Returns
        -------
        reconstructed_slab : pymatgen.core.surface.Slab
            Reconstructed slab with the desired parameters.

        """
        # Layers (containing sites) are removed from the bottom until
        # the desired slab_thickness is reached
        slab_resized = Geometry.remove_layers(slab, slab_thickness)
        
        # Necessary slab attributed to reconstruct the slab
        slab_attrs = ["species", "miller_index", "oriented_unit_cell", "shift", 
                      "scale_factor", "reorient_lattice", "reconstruction", 
                      "site_properties", "energy"]
        slab_params = attr_to_dict(slab_resized, slab_attrs)
        
        # To avoid issues with fractional coordinates, cartesian coordinates are used
        corrected_params = {'coords': slab_resized.cart_coords, 'coords_are_cartesian': True}
        slab_params.update(corrected_params)
        
        # Initial vacuum region is calculated
        initial_vacuum = Geometry.get_vacuum(slab_resized)
        
        # Lattice parameters are generated in order to be modified
        lat_attrs = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        lat_params = attr_to_dict(slab_resized.lattice, lat_attrs)
        
        proj_height = Geometry.get_proj_height(slab_resized)
        
        # 'c' parameter of the Lattice is modified to adjust vacuum to the desired thickness
        lat_params['c'] += (vacuum_thickness-initial_vacuum)*lat_params['c']/proj_height
        new_lat = Lattice.from_parameters(**lat_params)
        
        # Reconstructed slab is generated from the resized slab parameters and modified Lattice
        reconstructed_slab = center_slab(Slab(new_lat, **slab_params))
        return reconstructed_slab
        
    @staticmethod
    def get_surface_normal(struct):
        """
        Simple method to calculate the vector that is normal to the surface.
        Used to calculate the projected height of the structures. Surface
        normal is calculated in such a way that the surface is parallel to the
        first 2 lattice vectors.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            Main object in pymatgen to store structures.

        Returns
        -------
        numpy.ndarray
            3-dimensional numpy array that stores the surface normal.

        """
        latvec = struct.lattice.matrix
        normal = np.cross(latvec[0], latvec[1])
        return normal/np.linalg.norm(normal)
    
    @staticmethod
    def get_proj_height(struct):
        """
        Simple method to calculate the projected height of structures.
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
        surface_normal = Geometry.get_surface_normal(struct)
        return np.dot(struct.lattice.matrix[2], surface_normal)
    
    @staticmethod
    def get_layers(struct, tol = 0.1):
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
        
        # for correct layering, we need the distances between sites projected 
        # in the miller direction, this is necessary for nonorthogonal unit cells
        # proj_height = Geometry.get_proj_height(struct)
        
        # number of sites in the structure
        n = len(struct)
        frac_coords = struct.frac_coords
        
        # initiate a num_sites dimensional square distance matrix and populate it
        dist_matrix = np.zeros((n,n))
        for i, j in itertools.combinations(list(range(n)), 2):
            if i != j:
                cdist = frac_coords[i][2] - frac_coords[j][2]
                # cdist = abs(cdist - round(cdist)) * proj_height
                cdist = abs(cdist - round(cdist)) * struct.lattice.abc[2]
                dist_matrix[i, j] = cdist
                dist_matrix[j, i] = cdist
                
        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)
        
        # cluster the sites in the structure based on their c-coordinate and a given tolerance 
        clusters = fcluster(z, tol, criterion="distance")
        layers = defaultdict(list)
        for i,v in enumerate(clusters):
            layers[v].append(i)
            
        # for each layer, find sites that belong to it and assign an average c-value for the layer
        layers = {sum([struct.frac_coords[i][2] - np.floor(struct.frac_coords[i][2])
                       for i in v])/len(v):v for k,v in layers.items()}
        return layers
    
    @staticmethod
    def remove_layers(slab, num_layers, method = 'target', position = 'bottom'):
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
            Copy of the input Slab structure with layers removed
            and energy attribute modified.

        """
        layers = Geometry.get_layers(slab)
        if num_layers > len(layers):
            raise ValueError('Number of layers to remove/target can not exceed the number of layers \
                             in the given slab. Exiting..')
        c_coords = sorted(layers.keys())
        if method == "layers":
            to_remove = c_coords[:num_layers] if position == "bottom" else c_coords[num_layers:]
        elif method == "target":
            to_remove = c_coords[:len(c_coords)-num_layers] if position == "bottom" \
                else c_coords[len(c_coords)-num_layers:]
        indices_list = [layers[c_coord] for c_coord in to_remove]
        flat_list = [item for sublist in indices_list for item in sublist]
        slab_copy = slab.copy()
        slab_copy.remove_sites(flat_list)
        return slab_copy
    
def test_geometry(el, miller, slab_thick, vac_thick):
    # from triboflow.utils.database import NavigatorMP
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.structure import Structure
    # nav_mp = NavigatorMP()
    struct = Structure(Lattice.cubic(2.8), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    # struct, mp_id = nav_mp.get_low_energy_structure(el)
    bulk_conv = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    
    SG = SlabGenerator(initial_structure = bulk_conv,
                       miller_index = miller,
                       center_slab = True,
                       primitive = True,
                       in_unit_planes = True,
                       lll_reduce = False,
                       max_normal_search=max([abs(l) for l in miller]),
                       min_slab_size = slab_thick,
                       min_vacuum_size = vac_thick)
    
    slabs = SG.get_slabs(bonds=None, ftol=0.1, tol=0.1, symmetrize = False, repair=False)
    
    print('For {}-{} we have {} terminations'.format(struct.composition.reduced_formula, ''.join([str(i) for i in miller]), len(slabs)))
    print('Required slab thickness is {} layers and vacuum thickness is {} A'.format(slab_thick, vac_thick))
    print('resultant slabs have {} layers'.format([len(Geometry.get_layers(slab)) for slab in slabs]))
    print('resultant cells have {} A vacuum'.format([Geometry.get_vacuum(slab) for slab in slabs]))
    rec_slabs = [Geometry.reconstruct_slab(slab, slab_thick, vac_thick) for slab in slabs]
    print('Reconstructed slabs have {} layers'.format([len(Geometry.get_layers(slab)) for slab in rec_slabs]))
    print('Reconstructed cells have {} A vacuum'.format([Geometry.get_vacuum(slab) for slab in rec_slabs]))
    return slabs, rec_slabs