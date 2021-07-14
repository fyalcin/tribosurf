#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 16:39:05 2021

Class and methods that deal with matching interface structures.

The module contains:

    ** InterfaceMatcher **:
        Class to match two slabs to form an interface. The lattice search
        (_find_lattice_match method) is done by the implementation of the
        algorithm by Zur and McGill
        (Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084)
        as implemented in python in the MPInterfaces package
        (Computational Materials Science 122 (2016) 183–190;
        doi: 10.1016/j.commatsci.2016.05.020). However, in
        contrast to the MPInterfaces implementation the strain put on the
        lattices to get them to match is more flexible and achieved via a
        weighted average.
        The class contains the following modules, of which the last 4 alone
        usually are providing enough flexibility for the user.
            - __init__
            - __get_interface_dist
            - __flip_slab
            - __assign_top_bottom
            - __make_3d_lattice_from_2d_lattice
            - __get_supercell_matrix
            - _find_lattice_match
            - _get_matching_lattices
            - _get_matching_supercells
            - get_aligned_slabs
            - get_centered_slabs
            - get_interface
            - get_interface_distance

    Functions:
    - get_average_lattice
    - are_slabs_aligned

    Author: Michael Wolloch
            michael.wolloch@univie.ac.at

"""

import numpy as np
import warnings
from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.surface import center_slab, Slab
from mpinterfaces.transformations import get_matching_lattices

from triboflow.utils.structure_manipulation import (
    recenter_aligned_slabs,
    stack_aligned_slabs,
    clean_up_site_properties)
from triboflow.phys.shaper import Shaper


def get_average_lattice(latt1, latt2, weight1, weight2):
    """ Return a Lattice that is the weighted average of the input lattices.
    
    The weighted average of the lattice parameters and angles of two input
    lattices is returned. The orignal lattices are not changed. All lattices
    have to be pymatgen.core.lattice.Lattice objects.
    This is useful to get a single lattice for interface matching when two
    close lattices are found by another algorithm.
    

    Parameters
    ----------
    latt1 : pymatgen.core.lattice.Lattice
        First lattice for the averaging.
    latt2 : pymatgen.core.lattice.Lattice
        Second lattice for the averaging.
    weight1 : float
        Weight factor for the first lattice
    weight2 : float
        Weight factor for the first lattice

    Returns
    -------
    av_latt : pymatgen.core.lattice.Lattice
        Averaged lattice

    """
    a = np.average([latt1.a, latt2.a], weights=[weight1, weight2])
    b = np.average([latt1.b, latt2.b], weights=[weight1, weight2])
    c = np.average([latt1.c, latt2.c], weights=[weight1, weight2])
    alpha = np.average([latt1.alpha, latt2.alpha], weights=[weight1, weight2])
    beta = np.average([latt1.beta, latt2.beta], weights=[weight1, weight2])
    gamma = np.average([latt1.gamma, latt2.gamma], weights=[weight1, weight2])
    
    av_latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    return av_latt

def are_slabs_aligned(slab_1, slab_2, prec=12):
    """Check if two slabs have the same x and y lattice (z coordinate may differ)

    Parameters
    ----------
    slab_1 : pymatgen.core.surface.Slab
        First slab of two, which alignement is checked.
    slab_2 : pymatgen.core.surface.Slab
        Second slab of two, which alignement is checked.
    prec : int, optional
        Fixes the allowed deviation. The lattices are rounded to 'prec' digits.
        The default is 12.

    Returns
    -------
    bool
        Aligned or not

    """
    m1 = np.round(slab_1.lattice.matrix, prec)
    m2 = np.round(slab_2.lattice.matrix, prec)
    if (m1[0:2,:] == m2[0:2,:]).all():
        return True
    else:
        return False
 
def flip_slab(slab):
    """
    Flip the z coordinates of the input slab by multiplying all z-coords with -1.

    Parameters
    ----------
    slab : pymatgen.core.surface.Slab
       The input slab object flip

    Returns
    -------
    flipped_slab : pymatgen.core.surface.Slab
        The flipped slab

    """
    flip_matrix = np.array([[1., 0., 0.],
                            [0., 1., 0.],
                            [0., 0., -1.]])
    flipped_coords = np.dot(slab.cart_coords, flip_matrix)
        
    flipped_slab = Slab(lattice=slab.lattice,
                        species=slab.species,
                        coords=flipped_coords,
                        miller_index=slab.miller_index,
                        oriented_unit_cell=slab.oriented_unit_cell,
                        shift=slab.shift,
                        scale_factor=slab.scale_factor,
                        reconstruction=slab.reconstruction,
                        coords_are_cartesian=True,
                        site_properties=slab.site_properties)
    return center_slab(flipped_slab)

class InterfaceMatcher:
   
    def __init__(self,
                 slab_1,
                 slab_2,
                 strain_weight_1=1.0,
                 strain_weight_2=1.0,
                 max_area=100.0,
                 max_mismatch=0.01,
                 max_angle_diff=1.0,
                 r1r2_tol=0.02,
                 vacuum_thickness=15,
                 best_match='area',
                 interface_distance='auto'):
        """Initialize the InterfaceMatcher class
        
        If the strain weights sum up to zero (which is not allowed for the
        averaging process) they will be internally reset to 1 and a warning
        will be raised. Thus equal strain will be put on the slabs to form
        the interface. 
        

        Parameters
        ----------
        slab_1 : pymatgen.core.surface.Slab
            First slab for the interface matching
        slab_2 : pymatgen.core.surface.Slab
            Second slab for the interface matching
        strain_weight_1 : float, optional
            Weight for the lattice straining of the first slab. Combined with the
            'strain_weight_2' this determines which slab will be strained by
            how much to get unified lattice for the interface. The default is 1.0 
        strain_weight_2 : float, optional
            Weight for the lattice straining of the first slab. Combined with the
            'strain_weight_1' this determines which slab will be strained by
            how much to get unified lattice for the interface.  The default is 1.0
        max_area : float, optional
            Maximally allowed cross section area of the interface.
            The default is 100.0
        max_mismatch : float, optional
            Maximally allowed mismatch between the length of the lattice vectors
            in %. E.g. 0.01 is indicating a maximal 1% mismatch.
            The default is 0.01
        max_angle_diff : float, optional
            Maximally allowed angle mismatch between the cells in degrees.
            The default is 1.0
        r1r2_tol : float, optional
            Tolerance parameter related to the lattice search
            (abs(float(r1) * area1 - float(r2) * area2) / max_area <= r1r2_tol)
            The default is 0.02
        best_match : 'area' or 'mismatch', optional
            Determines if the algorithm returns the matched slabs with the
            smallest mismatch within "max_area" or the smalles area within the
            other tolerance parameters. The default is "area"
        interface_distance : int, float or str, optional
            Determines the distance between the matched slabs if an interface
            structure is returned. If the input can be transformed into a float,
            that float will be used as the distance. If not, the distance will
            be automatically computed as the average layer distance of the
            two slabs. The default is "auto"

        Returns
        -------
        None.

        """
        # handle inconsistent input
        if not sum((strain_weight_1, strain_weight_2)):
            warnings.warn('The sum of the weights is zero which is not allowed!\n'
                          'The strain weights will both be reset to 1 and strain '
                          'will be distributed evenly between the two materials.')
            strain_weight_1 = 1
            strain_weight_2 = 1
        if best_match not in ['area', 'mismatch']:
            warnings.warn('You have passed the "best_match" argument "{}", '
                         'which is not in ["area", "mismatch"]. It will be '
                         'set to "area" instead!'.format(best_match))
            best_match = 'area'
            
        self.match_params = {'max_area': max_area,
                             'max_mismatch': max_mismatch,
                             'max_angle_diff': max_angle_diff,
                             'r1r2_tol': r1r2_tol,
                             'best_match': best_match
                            }
        self.vacuum_thickness = vacuum_thickness
        self.aligned_top_slab = None
        self.aligned_bot_slab = None
        # Assign top and bottom slabs with strain weights.
        self.__assign_top_bottom(slab_1.copy(),
                                 slab_2.copy(),
                                 strain_weight_1,
                                 strain_weight_2)
        # Set interface distance
        self.__set_interface_dist(interface_distance)
        # Set the vacua for the centered slabs to ensure correct vacuum for the
        # interface
        self.__set_vacua()
        
    def __set_vacua(self):
        thickness_top = Shaper._get_proj_height(self.top_slab, 'slab')
        thickness_bot = Shaper._get_proj_height(self.bot_slab, 'slab')
        self.vacuum_top = thickness_bot + self.vacuum_thickness + self.inter_dist
        self.vacuum_bot = thickness_top + self.vacuum_thickness + self.inter_dist

    def __set_interface_dist(self, initial_distance):
        """
        Set the interface distance for the class instance.
        
        If the input can be transformed into a float,
        that float will be used as the distance. If not, the distance will
        be automatically computed as the average layer distance of the
        two slabs.

        Parameters
        ----------
        initial_distance : any type possible
            if transformable to a float, the float will be used, otherwise
            automatic computation is done.

        Returns
        -------
        None.

        """
        try:
            self.inter_dist = float(initial_distance)
        except:
            av_spacing_top = Shaper._get_average_layer_spacing(self.top_slab)
            av_spacing_bot = Shaper._get_average_layer_spacing(self.bot_slab)
            self.inter_dist = np.mean([av_spacing_top, av_spacing_bot])
                
    
    def __get_formula_and_miller(self, slab):
        """
        Return a string combination of a reduced formula and miller indices.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            input slab

        Returns
        -------
        str
            Reduced formula plus miller indices

        """
        f, _ = slab.composition.get_reduced_formula_and_factor()
        m = ''.join(str(s) for s in slab.miller_index)
        return f+m
    
    
    def __assign_top_bottom(self, slab_1, slab_2, weight_1, weight_2):
        """
        Assign top and bottom slab based on the formula and miller index.
        
        This is just for consistency since above and below are of course
        not meaningfull in a DFT context. Strain weights are assigned
        accordingly as well.

        Parameters
        ----------
        slab_1 : pymatgen.core.surface.Slab
            First slab.
        slab_2 : pymatgen.core.surface.Slab
            Second slab.
        weight_1 : float
            Weight for the strain distribution for slab 1.
        weight_2 : float
            Weight for the strain distribution for slab 2.

        Returns
        -------
        None.

        """
        n1 = self.__get_formula_and_miller(slab_1)
        n2 = self.__get_formula_and_miller(slab_2)
        opt1 = min(n1, n2)
        opt2 = max(n1, n2)
        
        if n1 == opt1 and n2 == opt2:
            self.top_slab = slab_1
            self.top_weight = weight_1
            self.bot_slab = slab_2
            self.bot_weight = weight_2
        else:
            self.top_slab = slab_2
            self.top_weight = weight_2
            self.bot_slab = slab_1
            self.bot_weight = weight_1
            
            
    def __make_3d_lattice_from_2d_lattice(self, slab, uv):
        """
        Takes a slab and adds its third lattice vector to the 2D lattice that is also passed.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            slab from which the third lattice vector is taken
        uv : [np.array, np.array]
            2D lattice

        Returns
        -------
        latt : pymatgen.core.lattice.Lattice
            Full 3D lattice

        """

        latt = Lattice(np.array([uv[0][:],
                                 uv[1][:],
                                 slab.lattice.matrix[2, :]
                                 ]
                                ))
        return latt

    def __get_supercell_matrix(self, slab, lattice):
        """
        Return a matrix that can be used to construct a supercell.

        Parameters
        ----------
        slab : pymatgen.core.surface.Slab
            Input slab for which a supercell should be constructed
        lattice : pymatgen.core.lattice.Lattice
            Lattice to which the input slab should be transformed.

        Returns
        -------
        sc_matrix : np.array
            Scaling matrix to construct a supercell.

        """
        _, __, sc_matrix = slab.lattice.find_mapping(lattice,
                                ltol=self.match_params['max_mismatch'],
                                atol=self.match_params['max_angle_diff'])
        sc_matrix[2] = np.array([0, 0, 1])
        return sc_matrix
    
    def _set_intended_vacuum(self, slab, vacuum):
        return Shaper._modify_vacuum(slab,
                                     vacuum,
                                     method='to_value',
                                     center=False)
        
    
    def _find_lattice_match(self):
        """
        Compute a matching lattice for heterogeneous interfaces.
        
        Compute the reduced matching lattice vectors for heterostructure
        interfaces as described in the paper by Zur and McGill:
        Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
        The function used is taken from the MPInterfaces package:
        Computational Materials Science 122 (2016) 183–190;
        doi: http://dx.doi.org/10.1016/j.commatsci.2016.05.020


        Returns
        -------
        uv_opt1 : [np.array, np.array]
            2D lattice for the first slab.
        uv_opt2 : [np.array, np.array]
            2D lattice for the second slab.

        """
        
        uv_opt1, uv_opt2 = get_matching_lattices(self.top_slab,
                                                 self.bot_slab,
                                                 **self.match_params)
        return uv_opt1, uv_opt2
    
    def _get_matching_lattices(self):
        """
        Find a matching lattice and return two Lattices if it is found.

        Returns
        -------
        None, None
            Return None twice if no match is found with the parameters.
        pymatgen.core.lattice.Lattice, pymatgen.core.lattice.Lattice
            Return two Lattice objects that are matched if a match is found.

        """
        uv_1, uv_2 = self._find_lattice_match()
        # Handle the possibility that no match is found
        if not uv_1 and not uv_2:
            return None, None
        top_latt = self.__make_3d_lattice_from_2d_lattice(self.top_slab, uv_1)
        bot_latt = self.__make_3d_lattice_from_2d_lattice(self.bot_slab, uv_2)
        
        return top_latt, bot_latt
    
    
    def _get_matching_supercells(self):
        """
        Return almost matched supercells of the input slabs

        Returns
        -------
        None, None
            Return None twice if no match is found with the parameters.
        pymatgen.core.surface.Slab, pymatgen.core.surface.Slab
            Two slabs that are almost matched, but still need to be strained.

        """
        top_latt, bot_latt = self._get_matching_lattices()
        # Handle the possibility that no match is found
        if not top_latt and not bot_latt:
            return None, None
        supercell_matrix_top = self.__get_supercell_matrix(self.top_slab, top_latt)
        supercell_matrix_bot = self.__get_supercell_matrix(self.bot_slab, bot_latt)
        supercell_top = self.top_slab.copy()
        supercell_bot = self.bot_slab.copy()
        supercell_top.make_supercell(supercell_matrix_top)
        supercell_bot.make_supercell(supercell_matrix_bot)
        return supercell_top, supercell_bot
            

    def get_aligned_slabs(self):
        """
        Get alinged slabs that are stretched according to the supplied weights.
        
        The first slab returned will be flipped horizontally so that the side
        facing the interface will be the one initially facing in the positive
        z direction.

        Returns
        -------
        None, None
            Return None twice if no match is found with the parameters.
        pymatgen.core.surface.Slab, pymatgen.core.surface.Slab
            Two slabs that are fully matched, e.g. they have the same lattice
            in x and y directions, z lattice vector may differ though!

        """
        if are_slabs_aligned(self.top_slab, self.bot_slab):
            print("\n  Slabs are already aligned!\n")
            flipped_slab = center_slab(flip_slab(self.top_slab))
            self.aligned_top_slab, self.aligned_bot_slab = flipped_slab, self.bot_slab
            return self.aligned_top_slab, self.aligned_bot_slab
        else:
            sc_top, sc_bot = self._get_matching_supercells()
            # Return None, None if no match was found for the given parameters.
            if not sc_top and not sc_bot:
                return None, None
            new_lattice = get_average_lattice(sc_top.lattice,
                                              sc_bot.lattice,
                                              self.top_weight,
                                              self.bot_weight)
            l_top = self.__make_3d_lattice_from_2d_lattice(sc_top, new_lattice.matrix)
            l_bot = self.__make_3d_lattice_from_2d_lattice(sc_bot, new_lattice.matrix)
                                                        
            sc_top.lattice = l_top
            sc_bot.lattice = l_bot
            
            flipped_slab = center_slab(flip_slab(sc_top))
            self.aligned_top_slab, self.aligned_bot_slab = flipped_slab, sc_bot
            return self.aligned_top_slab, self.aligned_bot_slab        
        
    def get_centered_slabs(self):
        """
        Return slabs that are already positioned to form and interface around z=0.
        
        
        The first slab returned will be flipped horizontally so that the side
        facing the interface will be the one initially facing in the positive
        z direction.

        Returns
        -------
        None, None
            Return None twice if no match is found with the parameters.
        pymatgen.core.surface.Slab, pymatgen.core.surface.Slab
            Two slabs that are fully matched, e.g. they have the same lattice
            in x and y directions, z lattice vector may differ though! They
            are positioned in such a way that the first slab will be half the
            interface distance above z=0, while the second one will be
            half that distance below z=0.

        """
        if self.aligned_top_slab and self.aligned_bot_slab:
            top_slab, bot_slab = self.aligned_top_slab, self.aligned_bot_slab
        else:
            top_slab, bot_slab = self.get_aligned_slabs()
        if not top_slab and not bot_slab:
                return None, None
        top_slab = self._set_intended_vacuum(top_slab, self.vacuum_top)
        bot_slab = self._set_intended_vacuum(bot_slab, self.vacuum_bot)
        tcs, bcs = recenter_aligned_slabs(top_slab, bot_slab, d=self.inter_dist)
        return tcs, bcs
        
    def get_interface(self):
        """
        Return a interface slab object containing two matched slabs.
        
        They are correctly orientated and have a defined distance to each other.
        The interface plane is z=0. The relative lateral position of the slabs
        are random though.

        Returns
        -------
        pymatgen.core.surface.Slab
            Matched interface structure of two slabs.

        """
        tcs, bcs = self.get_centered_slabs()
        if not tcs and not bcs:
                return None
        interface = stack_aligned_slabs(bcs, tcs)
        clean_interface = clean_up_site_properties(interface)
        
        return clean_interface
    
    def get_interface_distance(self):
        """
        Return the interface distance of the InterfaceMatcher class instance.

        Returns
        -------
        float
            Distance between the two slabs in Angstrom.

        """
        return self.inter_dist