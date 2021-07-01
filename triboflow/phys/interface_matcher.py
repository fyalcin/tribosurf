#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 16:39:05 2021

@author: mwo
"""
import numpy as np
from pymatgen.core.lattice import Lattice
from mpinterfaces.transformations import get_matching_lattices

class MatchInterface:
   
    def __init__(self,
                 slab_1,
                 slab_2,
                 strain_weight_1,
                 strain_weight_2,
                 max_area,
                 max_mismatch,
                 max_angle_diff,
                 r1r2_tol):
        self.slab_1 = slab_1.copy()
        self.slab_2 = slab_2.copy()
        self.w1 = strain_weight_1
        self.w2 = strain_weight_2
        self.match_constrains = {'max_area': max_area,
                                 'max_mismatch': max_mismatch,
                                 'max_angle_diff': max_angle_diff,
                                 'r1r2_tol': r1r2_tol
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
                                                 **self.match_constrains)
        return uv_opt1, uv_opt2
    
    def __make_3d_lattice_from_2d_lattice(self, structure, uv):
        latt = Lattice(np.array([uv[0][:],
                                 uv[1][:],
                                 structure.lattice.matrix[2, :]
                                 ]
                                ))
        return latt
    
    def __get_supercell_matrix(self, structure, lattice):
        _, __, scell = structure.lattice.find_mapping(lattice,
                                ltol=self.match_constrains['max_mismatch'],
                                atol=self.match_constrains['max_angle_diff'])
        scell[2] = np.array([0, 0, 1])
        return scell
    
    def _get_matching_lattices(self):
        uv_1, uv_2 = self._find_lattice_match()
        latt_1 = self.__make_3d_lattice_from_2d_lattice(self.slab_1, uv_1)
        latt_2 = self.__make_3d_lattice_from_2d_lattice(self.slab_2, uv_2)
        return latt_1, latt_2
    
    
    def _get_matching_supercells(self):
        latt_1, latt_2 = self._get_matching_lattices()
        supercell_matrix_1 = self.__get_supercell_matrix(self.slab_1, latt_1)
        supercell_matrix_2 = self.__get_supercell_matrix(self.slab_2, latt_2)
        supercell_1 = self.slab_1.copy()
        supercell_2 = self.slab_2.copy()
        supercell_1.make_supercell(supercell_matrix_1)
        supercell_2.make_supercell(supercell_matrix_2)
        return supercell_1, supercell_2
    
    
    @staticmethod
    def get_average_lattices(latt1, latt2, weight1, weight2):

        a = np.average([latt1.a, latt2.a], weights=[weight1, weight2])
        b = np.average([latt1.b, latt2.b], weights=[weight1, weight2])
        c = np.average([latt1.c, latt2.c], weights=[weight1, weight2])
        alpha = np.average([latt1.alpha, latt2.alpha], weights=[weight1, weight2])
        beta = np.average([latt1.beta, latt2.beta], weights=[weight1, weight2])
        gamma = np.average([latt1.gamma, latt2.gamma], weights=[weight1, weight2])
        
        av_latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        return av_latt        
        

    def get_aligned_slabs(self):
        ss1, ss2 = self._get_matching_supercells()
        new_lattice = self.get_averaged_lattice(ss1.lattice,
                                                ss2.lattice,
                                                self.w1,
                                                self.w2)
        l1 = self.__make_3d_lattice_from_2d_lattice(ss1, new_lattice.matrix)
        l2 = self.__make_3d_lattice_from_2d_lattice(ss2, new_lattice.matrix)
        ss1.lattice = l1
        ss2.lattice = l2
        return ss1, ss2