#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:40:07 2022

@author: wolloch
"""

import numpy as np
from itertools import combinations


from triboflow.phys.high_symmetry import (
    get_slab_hs, get_interface_hs, pbc_hspoints, fix_hs_dicts)
from triboflow.phys.interface_matcher import flip_slab

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.structure_matcher import StructureMatcher


class HighSymmetryAnalysis:
    
    def __init__(self,
                 interface,
                 in_frac_coordinates = True,
                 no_obtuse_hollow = True):
        
        self.interface = interface
        self.top_slab = flip_slab(interface.film)
        self.bot_slab = interface.substrate
        
        self.top_adsf = AdsorbateSiteFinder(self.top_slab)
        self.bot_adsf = AdsorbateSiteFinder(self.bot_slab)
        
        self.frac_coords = in_frac_coordinates
        self.no_optuse_hollow = no_obtuse_hollow

    def __get_adsorption_site(self, slab, unique): 
        """
        Return the adorption sites of a slab object, either only unique ones, or all.
        
        Parameters
        ----------
        slab : pymatgen.core.surface.Slab (pymatgen.core.structure.Structure)
            The slab of which you want to calculate the surface HS points.
            The input need to be a slab, i.e. an atomic structure non-periodic
            along the z-direction. However, type(slab) could be either a Slab
            object or a Structure object.
        
        unique : bool, optional
            If True, only unique non equivalent points are returned. If Flase, all
            replicas are listed as well. The default is True

        Returns
        -------
        sites : dict
            Dictionary with adsorption sites with the sites ordered by type.
            
        """

        if slab == self.top_slab:
            adsf = self.top_adsf
        else:
            adsf = self.bot_adsf
    
        if unique:
            # Extract the unique HS points for the given surface
            sites = adsf.find_adsorption_sites(distance=0, 
                                               symm_reduce=0.01, 
                                               near_reduce=0.01,
                                               no_obtuse_hollow=self.no_optuse_hollow)
        else:
            #Extract all the HS points for the given surface
            sites = adsf.find_adsorption_sites(distance=0, 
                                               symm_reduce=0, 
                                               near_reduce=0.01,
                                               no_obtuse_hollow=self.no_optuse_hollow)
        return sites
            
    def __return_frac_coords(self, hs_dict):
        frac_sites = {}
        for k, v in hs_dict.items():
            try:
                frac_sites[k] = np.dot(v, self.interface.lattice.inv_matrix)
            except:
                frac_sites[k] = v
        return frac_sites
    
    def __remove_c_coord(self, hs_dict):
        new_dict = {}
        for k, v in hs_dict.items():
            try:
                new_dict[k] = v[:,:2]
            except:
                new_dict[k] = v[:2]
        return new_dict
    
    def __separate_adsorption_sites(self, sites_dictionary):
        d = {}
        for k, v in sites_dictionary.items():
            if k != 'all':
                for i, coords in enumerate(v):
                    d[k+'_'+str(i+1)] = coords
        return d
    
    def __sort_all_sites(self, adsf, all_sites, unique_sites):
        new_dict = {}
        for label, unique_site in unique_sites.items():
            new_dict[label] = []
            for site in all_sites['all']:
                if len(adsf.symm_reduce([site, unique_site])) == 1:
                    new_dict[label].append(site)
            new_dict[label] = np.asarray(new_dict[label])
        return new_dict
    
    def _get_slab_hs_dicts(self):
        top_u = self.__separate_adsorption_sites(
                    self.__get_adsorption_site(self.top_slab, unique=True))
        if self.frac_coords:
            self.top_sites_unique = self.__remove_c_coord(
                self.__return_frac_coords(top_u))
        else:
            self.top_sites_unique = self.__remove_c_coord(top_u)
            
        bot_u = self.__separate_adsorption_sites(
            self.__get_adsorption_site(self.bot_slab, unique=True))
        if self.frac_coords:
            self.bot_sites_unique = self.__remove_c_coord(
                self.__return_frac_coords(bot_u))
        else:
            self.bot_sites_unique = self.__remove_c_coord(bot_u)
        
        top_all = self.__sort_all_sites(
            self.top_adsf,
            self.__get_adsorption_site(self.top_slab, unique=False),
            top_u)
        if self.frac_coords:
            self.top_sites_all = self.__remove_c_coord(
                self.__return_frac_coords(top_all))
        else:
            self.top_sites_all = self.__remove_c_coord(top_all)
            
        bot_all = self.__sort_all_sites(
            self.bot_adsf,
            self.__get_adsorption_site(self.bot_slab, unique=False),
            bot_u)
        if self.frac_coords:
            self.bot_sites_all = self.__remove_c_coord(
                self.__return_frac_coords(bot_all))
        else:
            self.bot_sites_all = self.remove_c_coord(bot_all)
        
    def __remove_equivalent_shifts(self, hs_shifts):
        unique_shifts = hs_shifts.copy()
        struct_match = StructureMatcher(ltol=0.01, stol=0.01, angle_tol=0.01,
                                        primitive_cell=False, scale=False)
        #auxiliary list since we cannot remove entries of the dictionary while
        #iterating over it.
        remove_these_shifts = []
        interface_1 = self.interface.copy()
        interface_2 = self.interface.copy()
        for shift_combos in combinations(hs_shifts.items(), 2):
            interface_1.in_plane_offset = shift_combos[0][1]
            interface_2.in_plane_offset = shift_combos[1][1]
            #if we find equivalent structures for different shifts, remove one of them
            if struct_match.fit(interface_1, interface_2):
                if shift_combos[1][0] not in remove_these_shifts:
                    remove_these_shifts.append(shift_combos[1][0])
        #finally remove the shifts
        print(remove_these_shifts)
        for shift in remove_these_shifts:
            unique_shifts.pop(shift)
        return unique_shifts
                    
        
    def _get_unique_shifts(self):
        unique_shifts = {}
        for kbot, vbot in self.bot_sites_unique.items():
            for ktop, vtop in self.top_sites_unique.items():
                # fold back all coordinates into the unit cell
                x_shift = clean_frac_coord(np.round(vbot - vtop, 12)[0])
                y_shift = clean_frac_coord(np.round(vbot - vtop, 12)[1])
                shift = np.asarray([x_shift, y_shift])
                # only add a new shift if the exact same one is not already present
                if not any((shift == x).all() for x in list(unique_shifts.values())):
                    unique_shifts[kbot+'-'+ktop] = shift
        
        return self.__remove_equivalent_shifts(unique_shifts)

def clean_frac_coord(coord):
    while coord < 0.0:
        coord += 1.0
    while coord > 1.0:
        coord -= 1.0
    if coord == 1.0 or coord == 0.0:
        return 0.0
    else:
        return np.round(coord, 12)

def old_symm(interface):
    top_slab = interface.film
    bot_slab = interface.substrate
    
    # Top slab needs to be flipped to find the high symmetry points at the
    # interface.
    flipped_top = flip_slab(top_slab)
    top_hsp_unique, top_hsp_all = get_slab_hs(flipped_top)

    bottom_hsp_unique, bottom_hsp_all = get_slab_hs(bot_slab)

    cell = interface.lattice.matrix

    hsp_unique = get_interface_hs(bottom_hsp_unique, top_hsp_unique, cell)
    hsp_all = get_interface_hs(bottom_hsp_all, top_hsp_all, cell)

    c_hsp_u, c_hsp_a = fix_hs_dicts(hsp_unique, hsp_all,
                                    top_slab, bot_slab)

    b_hsp_u = pbc_hspoints(bottom_hsp_unique, cell)
    b_hsp_a = pbc_hspoints(bottom_hsp_all, cell)
    t_hsp_u = pbc_hspoints(top_hsp_unique, cell)
    t_hsp_a = pbc_hspoints(top_hsp_all, cell)
    
    return {'c_hsp_u': c_hsp_u,
            'c_hsp_a': c_hsp_a,
            'b_hsp_u': b_hsp_u,
            'b_hsp_a': b_hsp_a,
            't_hsp_u': t_hsp_u,
            't_hsp_a': t_hsp_a}


def new_symm(interface):
    top_slab = interface.film
    bot_slab = interface.substrate
    
    # Top slab needs to be flipped to find the high symmetry points at the
    # interface.
    flipped_top = flip_slab(top_slab)