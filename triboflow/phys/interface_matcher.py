#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 16:39:05 2021

@author: mwo
"""
import numpy as np
import warnings
from pymatgen.core.lattice import Lattice
from mpinterfaces.transformations import get_matching_lattices


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

class MatchInterface:
   
    def __init__(self,
                 slab_1,
                 slab_2,
                 strain_weight_1=1.0,
                 strain_weight_2=1.0,
                 max_area=100.0,
                 max_mismatch=0.01,
                 max_angle_diff=1.0,
                 r1r2_tol=0.02,
                 best_match='area'):
        """Initialize the MatchInterface class
        
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
            
        self.slab_1 = slab_1.copy()
        self.slab_2 = slab_2.copy()
        self.w1 = strain_weight_1
        self.w2 = strain_weight_2
        self.match_params = {'max_area': max_area,
                                 'max_mismatch': max_mismatch,
                                 'max_angle_diff': max_angle_diff,
                                 'r1r2_tol': r1r2_tol,
                                 'best_match': best_match
                                 }

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
        
        uv_opt1, uv_opt2 = get_matching_lattices(self.slab_1,
                                                 self.slab_2,
                                                 **self.match_params)
        return uv_opt1, uv_opt2
    
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
        latt_1 = self.__make_3d_lattice_from_2d_lattice(self.slab_1, uv_1)
        latt_2 = self.__make_3d_lattice_from_2d_lattice(self.slab_2, uv_2)
        return latt_1, latt_2
    
    
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
        latt_1, latt_2 = self._get_matching_lattices()
        # Handle the possibility that no match is found
        if not latt_1 and not latt_2:
            return None, None
        supercell_matrix_1 = self.__get_supercell_matrix(self.slab_1, latt_1)
        supercell_matrix_2 = self.__get_supercell_matrix(self.slab_2, latt_2)
        supercell_1 = self.slab_1.copy()
        supercell_2 = self.slab_2.copy()
        supercell_1.make_supercell(supercell_matrix_1)
        supercell_2.make_supercell(supercell_matrix_2)
        return supercell_1, supercell_2
            

    def get_aligned_slabs(self):
        """
        Get alinged slabs that are stretched according to the supplied weights.

        Returns
        -------
        None, None
            Return None twice if no match is found with the parameters.
        pymatgen.core.surface.Slab, pymatgen.core.surface.Slab
            Two slabs that are fully matched, e.g. they have the same lattice
            in x and y directions, z lattice vector may differ though!

        """
        ss1, ss2 = self._get_matching_supercells()
        # Return None, None if no match was found for the given parameters.
        if not ss1 and not ss2:
            return None, None
        new_lattice = get_average_lattice(ss1.lattice,
                                          ss2.lattice,
                                          self.w1,
                                          self.w2)
        l1 = self.__make_3d_lattice_from_2d_lattice(ss1, new_lattice.matrix)
        l2 = self.__make_3d_lattice_from_2d_lattice(ss2, new_lattice.matrix)
        ss1.lattice = l1
        ss2.lattice = l2
        return ss1, ss2