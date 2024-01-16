#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:40:07 2022

@author: wolloch
"""
from copy import deepcopy
from typing import Union

import numpy as np
from monty.json import jsanitize
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.interface import Interface
from pymatgen.core.surface import Slab

from htflow_utils.shaper import Shaper
from triboflow.utils.structure_manipulation import flip_slab


class InterfaceSymmetryAnalyzer:
    """
    This class analyses an Interface with respect to its lateral shifts
    in order to efficiently construct a potential energy surface (PES),
    also called generalized stacking fault or gamma-surface.

    The procedure is as follows:
        1) Find the high symmetry adsorption sites (ontop, bridge, hollow) for
        both the top layer of the bottom part of the interface and the
        bottom layer for the top part of the interface.
        2) Sort these sites according the different types (e.g. ontop_1,
        bridge_2, ...) and distinguish between unique sites and equivalent
        replicas.
        3) Combine the ordered sites to possible lateral shifts between the
        top and bottom parts of the interface and label them
        (e.g. ontop_2-hollow_3).
        4) Separate the unique shifts from equivalent ones by using a
        pymatgen StructureMatcher to group the resulting interfaces.
        5) Return a dictionary with the results.

    After initializing the class with an Interface object, you can use the
    method get_high_symmetry_info() to get a dictionary with the following
    keys and values:

        'unique_shifts', -> Unique combinations of the interfaces high symmetry
                            points. Use this shifts to calculate the energies.
        'all_shifts',    -> All combinations of the interfaces high symmetry
                            points. These fill the whole unit cell. Copy
                            the calculated energies from the unique points
                            over to this get a smooth PES interpolation.
        'top_high_symm_points_unique',     -> unique high symmetry points
                                              for the flipped top slab.
        'top_high_symm_points_all',        -> all high symmetry points
                                              for the flipped top slab.
        'bottom_high_symm_points_unique',  -> unique high symmetry points
                                              for the bottom slab.
        'bottom_high_symm_points_all'      -> unique high symmetry points
                                              for the bottom slab.

    Furthermore, you can use the method get_all_high_symmetry_interfaces()
    to get a list of all interfaces with unique shifts and individually optimized
    interface gaps (if you not pass a specific interface distance).

    :param interface: A pymatgen Interface object.
    :type interface: pymatgen.core.interface.Interface

    :param in_cartesian_coordinates: Return the high symmetry points and interface shifts in cartesian rather than
    fractional coordinates. The default is False.
    :type in_cartesian_coordinates: bool, optional

    :param no_obtuse_hollow: Selects if you want to add obtuse hollows to the high symmetry points. The default for
    this in pymatgen's AdsorbateSiteFinder is "True", which means that no obtuse hollows are added! Be careful when
    adding these sites, since it will result most likely in extremely many unique shifts and sample the unit cell
    very densely with associated huge computational cost. The default is True.
    :type no_obtuse_hollow: bool, optional

    :param jsanitize_output: If true, will run monty sanitize on the final dictionary so that the numpy arrays will
    be converted to lists and can be saved in a mongoDB database.
    :type jsanitize_output: bool, optional

    :param ltol: Fractional length tolerance for the StructureMatcher. Keep this low, since matching structures
    should match exactly. The default is 0.01.
    :type ltol: float, optional

    :param stol: Site tolerance for the StructureMatcher. Defined as the fraction of the average free length per
    atom. Keep this low, since matching structures should match exactly. The default is 0.01.
    :type stol: float, optional

    :param angle_tol: Angle tolerance for the StructureMatcher in degrees. Keep this low, since matching structures
    should match exactly. The default is 0.01.
    :type angle_tol: float, optional
    """

    def __init__(
            self,
            interface: Interface,
            in_cartesian_coordinates: bool = False,
            no_obtuse_hollow: bool = True,
            jsanitize_output: bool = True,
            ltol: float = 0.01,
            stol: float = 0.01,
            angle_tol: float = 0.01,
    ):
        """
        Initializes the InterfaceSymmetryAnalyzer class.

        :param interface: A pymatgen Interface object.
        :type interface: pymatgen.core.interface.Interface

        :param in_cartesian_coordinates: Return the high symmetry points and interface shifts in cartesian rather
        than fractional coordinates. The default is False.
        :type in_cartesian_coordinates: bool, optional

        :param no_obtuse_hollow: Selects if you want to add obtuse hollows to the high symmetry points. The default
        for this in pymatgen's AdsorbateSiteFinder is "True", which means that no obtuse hollows are added! Be
        careful when adding these sites, since it will result most likely in extremely many unique shifts and sample
        the unit cell very densely with associated huge computational cost. The default is True.
        :type no_obtuse_hollow: bool, optional

        :param jsanitize_output: If true, will run monty sanitize on the final dictionary so that the numpy arrays
        will be converted to lists and can be saved in a mongoDB database.
        :type jsanitize_output: bool, optional

        :param ltol: Fractional length tolerance for the StructureMatcher. Keep this low, since matching structures
        should match exactly. The default is 0.01.
        :type ltol: float, optional

        :param stol: Site tolerance for the StructureMatcher. Defined as the fraction of the average free length per
        atom. Keep this low, since matching structures should match exactly. The default is 0.01.
        :type stol: float, optional

        :param angle_tol: Angle tolerance for the StructureMatcher in degrees. Keep this low, since matching
        structures should match exactly. The default is 0.01.
        :type angle_tol: float, optional
        """

        self.interface = interface
        self.cart_coords = in_cartesian_coordinates
        self.no_obtuse_hollow = no_obtuse_hollow
        self.jsanitize_output = jsanitize_output
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        self.high_symmetry_info = None

    def __convert_to_frac(self, hs_dict: dict, only_2d: bool = False) -> dict:
        """
        Convert cartesian coordinates in a high_symmetry dictionary to fractional.
        """
        if only_2d:
            m = self.interface.lattice.inv_matrix[:2, :2]
        else:
            m = self.interface.lattice.inv_matrix

        frac_sites = {}
        for k, v in hs_dict.items():
            try:
                frac_sites[k] = np.dot(v, m)
            except:
                frac_sites[k] = v
        return frac_sites

    def __convert_to_cart(self, hs_dict: dict, only_2d: bool = True) -> dict:
        """
        Convert fractional coordinates in a high_symmetry dictionary to cartesian.
        """
        if only_2d:
            m = self.interface.lattice.matrix[:2, :2]
        else:
            m = self.interface.lattice.matrix

        cart_sites = {}
        for k, v in hs_dict.items():
            try:
                cart_sites[k] = np.dot(v, m)
            except:
                cart_sites[k] = v
        return cart_sites

    @staticmethod
    def __remove_c_coord(hs_dict: dict) -> dict:
        """
        Remove the c coordinate for all high symmetry points to make 2D shifts.
        """
        new_dict = {}
        for k, v in hs_dict.items():
            try:
                new_dict[k] = v[:, :2]
            except:
                new_dict[k] = v[:2]
        return new_dict

    def __get_adsorption_site(self, slab: Slab, unique: bool) -> dict:
        """
        Return the adsorption sites of a slab object, either only unique ones, or all.

        :param slab: pymatgen.core.surface.Slab (pymatgen.core.structure.Structure)
            The slab of which you want to calculate the surface HS points.
            The input needs to be a slab, i.e. an atomic structure non-periodic
            along the z-direction. However, type(slab) could be either a Slab
            object or a Structure object.
        :type slab: pymatgen.core.surface.Slab (pymatgen.core.structure.Structure)

        :param unique: If True, only unique non-equivalent points are returned. If False, all
            replicas are listed as well.
        :type unique: bool, optional

        :return: Dictionary with adsorption sites with the sites ordered by type.
        :rtype: dict
        """

        if slab == self.top_slab:
            adsf = self.top_adsf
        else:
            adsf = self.bot_adsf

        if unique:
            # Extract the unique HS points for the given surface
            sites = adsf.find_adsorption_sites(
                distance=0,
                symm_reduce=0.01,
                near_reduce=0.01,
                no_obtuse_hollow=self.no_obtuse_hollow,
            )
        else:
            # Extract all the HS points for the given surface
            sites = adsf.find_adsorption_sites(
                distance=0,
                symm_reduce=0,
                near_reduce=0.01,
                no_obtuse_hollow=self.no_obtuse_hollow,
            )
        return sites

    @staticmethod
    def __separate_adsorption_sites(sites_dictionary: dict) -> dict:
        """
        Label unique adsorption sites uniquely

        Pymatgen's AdsorptionSiteFinder will group e.g. two distinctly different
        ontop sites together. This method separates them into ontop_1 and ontop_2.

        """
        d = {}
        for k, v in sites_dictionary.items():
            if k != "all":
                for i, coords in enumerate(v):
                    d[k + "_" + str(i + 1)] = coords
        return d

    @staticmethod
    def __sort_all_sites(adsf: AdsorbateSiteFinder, all_sites: dict, unique_sites: dict) -> dict:
        """
        Sort duplicate high symmetry points wrt the labels of the unique points.

        Uses an instance of pymatgens adsorption site finder to assign all
        high symmetry sites in a dictionary to one of the previously defined
        unique sites (e.g. ontop_1, bridge_2, ...)

        Parameters
        ----------
        adsf : pymatgen.analysis.adsorption.AdsorptionSiteFinder
            The AdsorptionSiteFinder object associated with the slab in question.
        all_sites : dict
            Dictionary previously constructed containing all (also symmetrically
            equivalent) high symmetry surface sites of the slab in question.
        unique_sites : Dict
            Dictionary with the unique high symmetry points already sorted.

        Returns
        -------
        new_dict : dict
            Dictionary with the same keys as "unique_sites" but containing all
            replica points originally in all_sites.

        """
        new_dict = {}
        for label, unique_site in unique_sites.items():
            new_dict[label] = []
            for site in all_sites["all"]:
                if len(adsf.symm_reduce([site, unique_site])) == 1:
                    new_dict[label].append(site)
            new_dict[label] = np.asarray(new_dict[label])
        return new_dict

    def __get_unique_hs_sites(self, slab: Slab) -> tuple[dict, dict]:
        """
        Return unique adsorption sites for a slab.

        Both full cartesian coordinates and fractional ones without the
        z-direction are returned.
        """
        unique = self.__separate_adsorption_sites(
            self.__get_adsorption_site(slab, unique=True)
        )
        unique_frac = self.__remove_c_coord(self.__convert_to_frac(unique))
        return unique, unique_frac

    def __get_all_hs_sites(self, slab: Slab, unique_sites: dict, adsf: AdsorbateSiteFinder) -> dict:
        """
        Return all adsorption sites for a slab in fractional 2D coordinates.
        """
        all_sites = self.__sort_all_sites(
            adsf, self.__get_adsorption_site(slab, unique=False), unique_sites
        )
        all_sites = self.__remove_c_coord(self.__convert_to_frac(all_sites))
        return all_sites

    def __get_slab_hs_dicts(self):
        """
        Analyse the adsorption sites of the substrate and the flipped film.

        In the process high symmetry dictionaries are constructed for both
        the top and the bottom slabs, for unique and replicated adsorption sites.
        Cartesian coordinates are converted to fractional ones. Z-coordinates
        are deleted.
        """

        top_u, self.top_sites_unique = self.__get_unique_hs_sites(
            self.top_slab
        )
        bot_u, self.bot_sites_unique = self.__get_unique_hs_sites(
            self.bot_slab
        )

        self.top_sites_all = self.__get_all_hs_sites(
            self.top_slab, top_u, self.top_adsf
        )
        self.bot_sites_all = self.__get_all_hs_sites(
            self.bot_slab, bot_u, self.bot_adsf
        )

    def __get_all_shifts(self):
        """
        Makes a dictionary containing all unique lateral shifts of an interface.

        They are grouped according to high-symmetry point combinations, e.g.
        ontop_1-bridge_2.
        """
        shifts = {}
        # all_shifts_list will be used to store all previously encountered shifts
        # since we do not want to duplicate shifts.
        all_shifts_list = []
        for kbot, vbot in self.bot_sites_all.items():
            for ktop, vtop in self.top_sites_all.items():
                # we do need to check the dimensionality of the arrays to loop
                # over shifts and not the x and y coordinates in the single shift.
                bottom_shifts = [vbot] if vbot.ndim == 1 else vbot
                top_shifts = [vtop] if vtop.ndim == 1 else vtop
                # shift_list holds all shifts associated with an HSP combination
                shift_list = []
                for bs in bottom_shifts:
                    for ts in top_shifts:
                        # The %1 is used to map all fractional coordinates back
                        # into the unit cell and get rid of -0.0 etc. We convert
                        # to a list to facilitate comparison with previous shifts.
                        shift = np.round((bs - ts) % 1, 5).tolist()
                        # only add a new shift if the exact same one is not already present
                        if shift not in all_shifts_list:
                            all_shifts_list.append(shift)
                            shift_list.append(shift)
                # convert to numpy array and assign to HSP combination.
                # do not add empty lists
                if shift_list:
                    shifts[kbot + "-" + ktop] = np.asarray(shift_list)
        self.all_shifts = shifts

    def __group_structures(self):
        """
        Group unique and all interfacial shifts using a StructureMatcher

        Different shifts of the interface can result in symmetrically equivalent
        structures of the interface. These can be found by pymatgen's
        StructureMatcher using the .group_structures method.
        It returns a list of lists with equivalent structures. We take the first
        one of each list and put them in self.unique_shifts, while the whole
        list is put in self.replica_shifts.
        """
        interface_list = []
        struct_match = StructureMatcher(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            primitive_cell=False,
            scale=False,
        )
        for name, shift_array in self.all_shifts.items():
            for shift in shift_array:
                # make an interface with the current shift and set a name
                intrfc = self.interface.copy()
                intrfc.in_plane_offset = shift
                intrfc.name = name
                interface_list.append(intrfc)
        # group the interfaces in interface_list
        grouped_interfaces = struct_match.group_structures(interface_list)
        # set up the output dictionaries and populate them with the correct
        # shifts which are just the .in_plane_offset attributes of the grouped
        # interface objects. Assign the group nr as keys, so the replicated
        # points can be matched to the unique ones.
        unique_shifts = {}
        group_assignment = {}
        replic_shifts = {}
        for i, intrfc_group in enumerate(grouped_interfaces):
            unique_shifts["group_" + str(i + 1)] = intrfc_group[
                0
            ].in_plane_offset
            group_assignment["group_" + str(i + 1)] = intrfc_group[0].name
            shift_list = []
            for intrfc in intrfc_group:
                shift_list.append(intrfc.in_plane_offset)
            replic_shifts["group_" + str(i + 1)] = np.unique(
                np.asarray(shift_list), axis=0
            )

        self.unique_shifts = unique_shifts
        self.replica_shifts = replic_shifts
        self.group_assignment = group_assignment

    def __set_parameters(self, interface: Interface):
        """
        set a couple of parameters for the class instance depended on the
        interface passed to the __call__ method.
        """
        self.interface = interface
        self.top_slab = flip_slab(interface.film)
        self.bot_slab = interface.substrate

        self.top_adsf = AdsorbateSiteFinder(self.top_slab)
        self.bot_adsf = AdsorbateSiteFinder(self.bot_slab)

    def __check_cartesian_and_jsanitize(self, out_dict: dict) -> dict:
        """
        Return the output and possibly switch to cartesian coordinates and jsanitize.
        """
        if self.cart_coords:
            new_out = {}
            for k, v in out_dict.items():
                new_out[k] = self.__convert_to_cart(v)
        else:
            new_out = out_dict

        if self.jsanitize_output:
            return jsanitize(new_out)
        else:
            return new_out

    def get_high_symmetry_info(self):
        """
        Return information about the interfaces high-symmetry points and shifts

        Used to generate a potential energy surface, or generalized stacking
        fault surface, for the input interface. Information about the lateral
        shifts (in_plane_offset) and the high symmetry points for the substrate
        and film (flipped) are also returned.

        Returns
        -------
        dict
            A dictionary with the following entries:
                'unique_shifts', -> Unique combinations of the interfaces high
                                    symmetry points. Use this shifts to
                                    calculate the energies.
                'all_shifts',    -> All combinations of the interfaces high
                                    symmetry points. These fill the whole unit
                                    cell. Copy the calculated energies from the
                                    unique points over to this get a smooth PES
                                    interpolation.
                'top_high_symm_points_unique',     -> unique high symmetry points
                                                      for the flipped top slab.
                'top_high_symm_points_all',        -> all high symmetry points
                                                      for the flipped top slab.
                'bottom_high_symm_points_unique',  -> unique high symmetry points
                                                      for the bottom slab.
                'bottom_high_symm_points_all'      -> unique high symmetry points
                                                      for the bottom slab.

        """
        self.__set_parameters(self.interface)
        self.__get_slab_hs_dicts()
        self.__get_all_shifts()
        self.__group_structures()
        out_dict = {
            "unique_shifts": self.unique_shifts,
            "all_shifts": self.replica_shifts,
            "top_high_symm_points_unique": self.top_sites_unique,
            "top_high_symm_points_all": self.top_sites_all,
            "bottom_high_symm_points_unique": self.bot_sites_unique,
            "bottom_high_symm_points_all": self.bot_sites_all,
            "group_assignments": self.group_assignment,
        }
        self.high_symmetry_info = out_dict

        return self.__check_cartesian_and_jsanitize(out_dict)

    def get_all_high_symmetry_interfaces(
            self,
            interface_distance: Union[str, float] = None,
            distance_boost: float = 0.05,
            bond_dist_delta: float = 0.05
    ) -> list[Interface]:
        """
        Return a list of interfaces with high-symmetry lateral shifts.

        To fit a PES to the interface, it is necessary to have a set of
        interfaces with different lateral shifts. This method returns a list
        of interfaces with different shifts.
        At the same time the interface distance is increased until the
        minimum bond distance is larger than the minimum bond distance of the
        aligned slabs minus the bond_dist_delta parameter. This is to ensure
        that the interface is not too close to the slabs, which would lead to
        unphysical bond distances.

        Parameters
        ----------
        interface_distance: Union[str, float], optional
            If set to "auto" or None, the interface distance is increased until
            the minimum bond distance is larger than the minimum bond distance
            of the aligned slabs minus the bond_dist_delta parameter. If set to
            a float, the interface distance is set to that value.
        distance_boost : float, optional
            Boost the interface distance at each step by this amount.
            The default is 0.05.
        bond_dist_delta : float, optional
            Minimum bond distance of the interface is set to the minimum bond
            distance of the aligned slabs minus this value.
            The default is 0.05.

        Returns
        -------
        list
            List of Interface objects.

        """
        initial_interface = self.interface
        # set interface distance extremely low initially, since we want to
        # increase it until the minimum bond distance is equal to the
        # minimum bond distance of the aligned slabs minus the bond_dist_delta
        initial_interface.gap = 0.1  # in Angstrom

        bonds_top = Shaper.get_all_bonds(
            struct=initial_interface.film,
            r=min(initial_interface.film.lattice.abc),
        )
        bonds_bot = Shaper.get_all_bonds(
            struct=initial_interface.substrate,
            r=min(initial_interface.substrate.lattice.abc),
        )
        bonds = bonds_top + bonds_bot
        min_bond = min(bonds, key=lambda x: x[2])[2]

        if self.high_symmetry_info is None:
            hsp_dict = self.get_high_symmetry_info()
        else:
            hsp_dict = self.high_symmetry_info

        interfaces = []
        for group_name, shift in hsp_dict["unique_shifts"].items():
            interface = deepcopy(initial_interface)
            interface.in_plane_offset = shift
            # only update interface distance if it is set to auto
            # or if it is None. If a float is passed, use that.
            if interface_distance == "auto" or interface_distance is None:
                interface_bonds = Shaper.get_all_bonds(
                    struct=interface, r=min(interface.lattice.abc)
                )
                min_interface_bond = min(interface_bonds, key=lambda x: x[2])[2]

                while min_interface_bond <= min_bond - bond_dist_delta:
                    interface.gap = interface.gap + distance_boost
                    interface_bonds = Shaper.get_all_bonds(
                        struct=interface, r=min(interface.lattice.abc)
                    )
                    min_interface_bond = min(interface_bonds, key=lambda x: x[2])[
                        2
                    ]
            else:
                interface.gap = interface_distance

            interface.interface_properties.update(
                {
                    "high_symmetry_info": {
                        "group_name": group_name,
                        "has_shift": shift,
                    }
                }
            )
            interfaces.append(interface)
        return interfaces
