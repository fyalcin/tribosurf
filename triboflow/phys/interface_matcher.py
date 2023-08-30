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
        as implemented in python in the pymatgen. However, in
        contrast to the pymatgen implementation the strain put on the
        lattices to get them to match is more flexible and achieved via a
        weighted average.
        It is also possible to get a list of interfaces with different
        lateral shifts, which is useful for fitting a PES to the interface.
        The class contains the following modules, of which the last 5 alone
        usually are providing enough flexibility for the user.
            - __init__
            - __get_interface_dist
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
    - get_consolidated_comp_params
    - get_average_lattice
    - are_slabs_aligned

    Author: Michael Wolloch
            michael.wolloch@univie.ac.at

"""

import numpy as np
import warnings
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.core.interface import Interface
from pymatgen.core.lattice import Lattice

from hitmen_utils.db_tools import VaspDB
from hitmen_utils.shaper import Shaper
from triboflow.utils.structure_manipulation import (
    recenter_aligned_slabs,
    clean_up_site_properties,
)


def get_consolidated_comp_params(mpid1, mpid2, bulk_coll, db_file, high_level):
    nav = VaspDB(db_file=db_file, high_level=high_level)
    bulk_1 = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid1})
    bulk_2 = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid2})

    comp_params_1 = bulk_1["comp_parameters"]
    comp_params_2 = bulk_2["comp_parameters"]

    encut_1 = comp_params_1["encut"]
    encut_2 = comp_params_2["encut"]
    encut_inter = max(encut_1, encut_2)
    k_dens_1 = comp_params_1["k_dens"]
    k_dens_2 = comp_params_2["k_dens"]
    k_dens_inter = max(k_dens_1, k_dens_2)
    metal_1 = comp_params_1["is_metal"]
    metal_2 = comp_params_2["is_metal"]
    metal_inter = any((metal_1, metal_2))

    # take one of the full comp_params dict and update the values
    comp_params_inter = comp_params_1
    comp_params_inter["encut"] = encut_inter
    comp_params_inter["k_dens"] = k_dens_inter
    comp_params_inter["is_metal"] = metal_inter

    return comp_params_inter


def get_average_lattice(latt1, latt2, weight1, weight2):
    """Return a Lattice that is the weighted average of the input lattices.

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
        First slab of two, which alignment is checked.
    slab_2 : pymatgen.core.surface.Slab
        Second slab of two, which alignment is checked.
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
    if (m1[0:2, :] == m2[0:2, :]).all():
        return True
    else:
        return False


class InterfaceMatcher:
    def __init__(
        self,
        slab_1,
        slab_2,
        strain_weight_1=1.0,
        strain_weight_2=1.0,
        max_area_ratio_tol=0.09,
        max_area=100,
        max_length_tol=0.03,
        max_angle_tol=0.01,
        bidirectional_search=True,
        vacuum=15.0,
        interface_distance="auto",
        interface_distance_addon=0.0,
        pressure=None,
    ):
        """Initialize the InterfaceMatcher class

        If the strain weights sum up to zero (which is not allowed for the
        averaging process), or one or both of them are negative, they will be
        internally reset to 1 and a warning will be raised. In that case,
        equal strain will be put on the slabs to form the interface. If
        the default of "auto" is used for "interface_distance", the average
        layer distance of both slabs is used for the interface spacing.
        By default, a small cell is prioritized in the lattice search, but
        "best_match" can be set to "mismatch" to find the smallest possible
        mismatch while staying below "max_area".


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
        max_area_ratio_tol : float, optional
            Tolerance parameter related to the lattice search
            Lattice vectors are considered for the search if the areas do not
            differ more than max_area_ratio_tol
            The default is 0.09
        max_area : float, optional
            Maximally allowed cross-sectional area of the interface.
            The default is 100.0
        max_length_tol : float, optional
            Maximally allowed mismatch between the length of the lattice vectors
            in %. E.g. 0.01 is indicating a maximal 1% mismatch.
            The default is 0.03
        max_angle_tol : float, optional
            Maximally allowed angle mismatch between the cells in %.
            E.g. 0.01 is indicating a maximal 1% mismatch.
            The default is 0.01
        bidirectional_search : bool, optional
            Select if vectors of substrate and film are exchangeable during the
            lattice search. The default is True
        vacuum : float, optional
            Thickness of the vacuum layer for the final interface. The default
            is 15.0
        interface_distance : int, float or str, optional
            Determines the distance between the matched slabs if an interface
            structure is returned. If the input can be transformed into a float,
            that float will be used as the distance. If not, the distance will
            be automatically computed as the average layer distance of the
            two slabs. The default is "auto"
        interface_distance_addon : float, optional
            This will be added to the automatic interface distance if
            interface_distance is not set explicitly to a float or int.
            The default is 0.0
        pressure : float, optional
            External pressure to be applied to the interface structure in GPa.
            For now this is doing nothing. However, one could use this to
            reduce the automatic interface distance if it is not explicitly
            set. The default is None.

        Returns
        -------
        None.

        """
        # handle inconsistent input
        self.strain = None
        self.interface = None
        (
            weight_1,
            weight_2,
        ) = self.__check_weights(strain_weight_1, strain_weight_2)

        self.match_params = {
            "max_area": max_area,
            "max_length_tol": max_length_tol,
            "max_angle_tol": max_angle_tol,
            "max_area_ratio_tol": max_area_ratio_tol,
            "bidirectional": bidirectional_search,
        }
        self.vacuum_thickness = vacuum
        self.aligned_top_slab = None
        self.aligned_bot_slab = None
        # Assign top and bottom slabs with strain weights.
        self.__assign_top_bottom(slab_1.copy(), slab_2.copy(), weight_1, weight_2)
        # Set interface distance
        self.__set_interface_dist(interface_distance, interface_distance_addon)
        # Set the vacua for the centered slabs to ensure correct vacuum for the
        # interface
        self.__set_vacua()

    def __check_weights(self, strain_weight_1, strain_weight_2):
        """
        Check if the supplied strain weights are making sense and correct them if not.


        Parameters
        ----------
        strain_weight_1 : float
            Weight 1
        strain_weight_2 : float
            weight 2

        Returns
        -------
        float
            Weight 1
        float
            Weight 2

        """
        if not sum((strain_weight_1, strain_weight_2)):
            warnings.warn(
                "The sum of the weights is zero which is not allowed!\n"
                "The strain weights will both be reset to 1 and strain "
                "will be distributed evenly between the two materials."
            )
            return 1.0, 1.0
        elif strain_weight_1 < 0 or strain_weight_2 < 0:
            warnings.warn(
                "One or both of the weights is negative which may "
                "lead to unforeseen effects and is not allowed!\n"
                "The strain weights will both be reset to 1 and strain "
                "will be distributed evenly between the two materials."
            )
            return 1.0, 1.0
        else:
            return strain_weight_1, strain_weight_2

    def __set_vacua(self):
        """
        Set the vacuum thickness for the top and bottom slabs.

        Returns
        -------
        None.

        """
        thickness_top = Shaper.get_proj_height(self.top_slab, "slab")
        thickness_bot = Shaper.get_proj_height(self.bot_slab, "slab")
        self.vacuum_top = thickness_bot + self.vacuum_thickness + self.inter_dist
        self.vacuum_bot = thickness_top + self.vacuum_thickness + self.inter_dist

    def __set_interface_dist(self, initial_distance, distance_boost):
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
        distance_boost : float or int
            added to interface distance if it is automatically determined.

        Returns
        -------
        None.

        """
        try:
            self.inter_dist = float(initial_distance)
        except ValueError:
            av_spacing_top = Shaper.get_average_layer_spacing(self.top_slab)
            av_spacing_bot = Shaper.get_average_layer_spacing(self.bot_slab)
            self.inter_dist = np.mean([av_spacing_top, av_spacing_bot]) + distance_boost

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
        m = "".join(str(s) for s in slab.miller_index)
        return f + m

    def __assign_top_bottom(self, slab_1, slab_2, weight_1, weight_2):
        """
        Assign top and bottom slab based on the formula and miller index.

        Make sure that the c-axis are orthogonal to the ab-plane.
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
            self.top_slab = slab_1.get_orthogonal_c_slab()
            self.top_weight = weight_1
            self.top_miller = slab_1.miller_index
            self.bot_slab = slab_2.get_orthogonal_c_slab()
            self.bot_weight = weight_2
            self.bot_miller = slab_2.miller_index
        else:
            self.top_slab = slab_2.get_orthogonal_c_slab()
            self.top_weight = weight_2
            self.top_miller = slab_2.miller_index
            self.bot_slab = slab_1.get_orthogonal_c_slab()
            self.bot_weight = weight_1
            self.bot_miller = slab_1.miller_index

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

        latt = Lattice(np.array([uv[0][:], uv[1][:], slab.lattice.matrix[2, :]]))
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
        _, __, sc_matrix = slab.lattice.find_mapping(
            lattice, ltol=self.match_params["max_length_tol"], atol=1
        )
        sc_matrix[2] = np.array([0, 0, 1])
        return sc_matrix

    def _set_intended_vacuum(self, slab, vacuum):
        return Shaper.modify_vacuum(slab, vacuum, method="to_value", center=False)

    def _find_lattice_match(self):
        """
        Compute a matching lattice for heterogeneous interfaces.

        Compute the reduced matching lattice vectors for heterostructure
        interfaces as described in the paper by Zur and McGill:
        Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
        Implementation by pymatgen.


        Returns
        -------
        matched_top_vec : [np.array, np.array]
            2D lattice for the first slab.
        matched_bot_vec : [np.array, np.array]
            2D lattice for the second slab.

        """
        zsl_gen = ZSLGenerator(**self.match_params)
        try:
            zsl_match = next(
                zsl_gen(
                    film_vectors=self.top_slab.lattice.matrix[0:2],
                    substrate_vectors=self.bot_slab.lattice.matrix[0:2],
                    lowest=True,
                )
            )
            matched_top_vec = zsl_match.film_sl_vectors
            matched_bot_vec = zsl_match.substrate_sl_vectors
            self.zsl_match = zsl_match
            return matched_top_vec, matched_bot_vec
        except:
            return None, None

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
        top_vec, bot_vec = self._find_lattice_match()
        # Try to make 3D lattices out of the 2D lattices, but take care of
        # the possibility that no match is found and both top_vec and bot_vec
        # are None.
        try:
            top_latt = self.__make_3d_lattice_from_2d_lattice(self.top_slab, top_vec)
            bot_latt = self.__make_3d_lattice_from_2d_lattice(self.bot_slab, bot_vec)
        except TypeError:
            return None, None

        self.unstrained_top_lattice = top_latt
        self.unstrained_bot_lattice = bot_latt

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
            self.aligned_top_slab, self.aligned_bot_slab = (
                self.top_slab,
                self.bot_slab,
            )
        else:
            sc_top, sc_bot = self._get_matching_supercells()
            # Return None, None if no match was found for the given parameters.
            if not sc_top and not sc_bot:
                return None, None
            new_lattice = get_average_lattice(
                sc_top.lattice,
                sc_bot.lattice,
                self.top_weight,
                self.bot_weight,
            )
            l_top = self.__make_3d_lattice_from_2d_lattice(sc_top, new_lattice.matrix)
            l_bot = self.__make_3d_lattice_from_2d_lattice(sc_bot, new_lattice.matrix)

            sc_top.lattice = l_top
            sc_bot.lattice = l_bot
            self.aligned_top_slab, self.aligned_bot_slab = sc_top, sc_bot

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
        Return an Interface object containing two matched slabs.

        They are correctly orientated and have a defined distance to each other.
        The relative lateral position of the slabs are random though.

        Returns
        -------
        pymatgen.core.interface.Interface
            Matched interface structure of two slabs.

        """
        tcs, bcs = self.get_centered_slabs()
        if not tcs and not bcs:
            return None

        # Note that the from_slab method of the Inteface object flips the film over!
        # UPDATE ON 27.06.22: After noticing some discrepancies with the high symmetry points
        # we decided to leave the film as is, since it gets flipped once more by the Interface
        # class. This seems to fix the issues with the high symmetry points.
        self.interface = Interface.from_slabs(
            substrate_slab=bcs,
            film_slab=tcs,
            gap=self.inter_dist,
            vacuum_over_film=self.vacuum_thickness,
        )

        self.interface = clean_up_site_properties(self.interface)

        self.interface.interface_properties = {
            "area": tcs.surface_area,
            "strain": self.get_strain(),
            "film_miller": self.top_miller,
            "substrate_miller": self.bot_miller,
            "strain_weights": {
                "film": self.top_weight,
                "substrate": self.bot_weight,
            },
        }

        return self.interface

    def get_interface_distance(self):
        """
        Return the interface distance of the InterfaceMatcher class instance.

        Returns
        -------
        float
            Distance between the two slabs in Angstrom.

        """
        return self.inter_dist

    def get_strain(self):
        """
        Calculate the strains on the aligned top and bottom slabs by comparing
        the unstrained lattices of the constituent slabs and final lattice
        of the interface.

        Returns
        -------
        dict
            Dictionary with keys 'top' and 'bot' with the percentages of
            strain in the two directions of the lattice vectors.

        """
        if hasattr(self, "strain"):
            return self.strain

        if are_slabs_aligned(self.top_slab, self.bot_slab):
            self.strain = {"top": [0, 0], "bot": [0, 0]}
            return self.strain

        if not hasattr(self, "interface"):
            self.get_interface()

        u_top_latt = self.unstrained_top_lattice
        u_bot_latt = self.unstrained_bot_lattice
        inter_latt = self.interface.lattice

        top_strain_ab = [
            100 * (inter_latt.abc[i] - u_top_latt.abc[i]) / u_top_latt.abc[i]
            for i in range(2)
        ]
        bot_strain_ab = [
            100 * (inter_latt.abc[i] - u_bot_latt.abc[i]) / u_bot_latt.abc[i]
            for i in range(2)
        ]

        top_strain_gamma = inter_latt.gamma - u_top_latt.gamma

        bot_strain_gamma = inter_latt.gamma - u_bot_latt.gamma

        self.strain = {
            "top": {
                "a": top_strain_ab[0],
                "b": top_strain_ab[1],
                "gamma": top_strain_gamma,
            },
            "bot": {
                "a": bot_strain_ab[0],
                "b": bot_strain_ab[1],
                "gamma": bot_strain_gamma,
            },
        }

        return self.strain
