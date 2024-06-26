from typing import Union

import numpy as np
from pymatgen.core.interface import Interface
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, Slab, center_slab
from pymatgen.ext.matproj import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from htflow_utils.db_tools import VaspDB
from htflow_utils.misc_tools import transfer_average_magmoms
from triboflow.utils.mp_connection import MPConnection


def get_interface_distance(structure: Interface) -> float | None:
    z_max_substrate = -10000
    z_min_film = 10000
    for s in structure.sites:
        if s.properties["interface_label"] == "substrate":
            if s.coords[-1] > z_max_substrate:
                z_max_substrate = s.coords[-1]
        elif s.properties["interface_label"] == "film":
            if s.coords[-1] < z_min_film:
                z_min_film = s.coords[-1]
        else:
            return None
    return z_min_film - z_max_substrate


def get_conv_bulk_from_mpid(mpid: str,
                            coll: str,
                            db_file: str = "auto",
                            high_level: Union[str, bool] = "auto") -> Structure:
    """
    Generates the conventional standard bulk structure for a material given by
    its MaterialsProject ID.

    Parameters
    ----------
    mpid : str
        Unique MaterialsProject ID describing the structure.
    coll : str
        Collection in which to search for the primitive bulk structure
        in the database.
    db_file : str, optional
        Full path of the db.json file. If 'auto', looks it up from the
        FW_CONFIG_FILE environment variable. The default is 'auto'.
    high_level : str, bool optional
        High level database in which to search for the primitive bulk
        structure. If True, looks it up from the db.json file.
        The default is True.

    Raises
    ------
    ValueError
        Either the fully optimized primitive bulk structure, or the one directly
        loaded from MP should already be in the database.

    Returns
    -------
    bulk_conv : pymatgen.core.structure.Structure
        Conventional standard bulk structure for the material.

    """
    db = VaspDB(db_file, high_level)
    bulk_dict = db.find_data(coll, {"mpid": mpid})
    if bulk_dict is None:
        with MPRester() as mpr:
            bulk_conv = mpr.get_structure_by_material_id(
                material_id=mpid, conventional_unit_cell=True
            )
        return bulk_conv
        # raise ValueError(f'Bulk entry for {mpid} not found in {high_level}')

    bulk_struct = bulk_dict.get("structure_equiVol")
    if bulk_struct is None:
        print(
            f"Optimized bulk structure not found for {mpid}, falling back to"
            "MaterialsProject imported structure."
        )
        bulk_struct = bulk_dict.get("structure_fromMP")
        if bulk_struct is None:
            raise ValueError(
                f"Bulk structure not found in the bulk entry for {mpid}"
            )
    bulk_struct = Structure.from_dict(bulk_struct)
    bulk_conv = SpacegroupAnalyzer(
        bulk_struct
    ).get_conventional_standard_structure(keep_site_properties=True)
    return bulk_conv


def slab_from_structure(miller: list | tuple,
                        structure: Structure) -> Slab:
    """Returns a pymatgen.core.surface.Slab from a pymatgen structure.

    Parameters
    ----------
    miller : list of int
        Miller indices given as a list of integers
    structure : pymatgen.core.structure.Structure
        Structure to be converted into a Slab

    Returns
    -------
    pymatgen.core.surface.Slab
        The input structure converted to a Slab

    """
    return Slab(
        lattice=structure.lattice,
        species=structure.species_and_occu,
        coords=structure.frac_coords,
        miller_index=miller,
        oriented_unit_cell=structure,
        shift=0,
        scale_factor=np.eye(3, dtype=int),
        site_properties=structure.site_properties,
    )


def clean_up_site_properties(structure: Structure) -> Structure:
    """
    Cleans up site_properties of structures that contain NoneTypes.

    If an interface is created from two different structures, it is possible
    that some site properties like magmom are not set for both structures.
    This can lead later to problems since they are replaced by None.
    This function replaces NoneTypes with 0.0 for magmom and deletes all other
    site_properties if None entries are found in it.


    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Input structure

    Returns
    -------
    struct : pymatgen.core.structure.Structure
        Output structure

    """
    struct = structure.copy()
    for key in struct.site_properties.keys():
        if key == "magmom":
            new_magmom = []
            for m in struct.site_properties[key]:
                if m is None:
                    new_magmom.append(0.0)
                else:
                    new_magmom.append(m)
            struct.add_site_property("magmom", new_magmom)
        else:
            if any(struct.site_properties[key]) is None:
                struct.remove_site_property(key)
    return struct


def stack_aligned_slabs(bottom_slab: Slab,
                        top_slab: Slab,
                        top_shift: tuple = (0, 0, 0)) -> Structure:
    """
    Combine slabs that are centered around 0 into a single structure.

    Optionally shift the top slab by a vector of cartesian coordinates.

    Parameters
    ----------
    bottom_slab : pymatgen.core.structure.Structure or pymatgen.core.surface.Slab
        Bottom slab.
    top_slab : pymatgen.core.structure.Structure or pymatgen.core.surface.Slab
        Top slab.
    top_shift : list of 3 floats, optional
        Vector of cartesian coordinates with which to shift the top slab.
        The default is [0,0,0].

    Returns
    -------
    interface : pymatgen.core.structure.Structure or pymatgen.core.surface.Slab
                depending on type of bottom_slab
        An interface structure of two slabs with an optional shift of the top
        slab.

    """
    interface = bottom_slab.copy()
    t_copy = top_slab.copy()

    t_copy.translate_sites(
        indices=range(len(t_copy.sites)),
        vector=top_shift,
        frac_coords=False,
        to_unit_cell=False,
    )

    for s in t_copy.sites:
        new_site = PeriodicSite(
            lattice=interface.lattice,
            coords=s.frac_coords,
            coords_are_cartesian=False,
            species=s.species,
            properties=s.properties,
        )
        interface.sites.append(new_site)

    return interface


def recenter_aligned_slabs(top_slab: Slab,
                           bottom_slab: Slab,
                           d: float = 2.5):
    """
    Center two slabs around z=0 and give them the distance d.

    Parameters
    ----------
    top_slab : pymatgen.core.structure.Structure
        The slab that should be on top.
    bottom_slab : pymatgen.core.structure.Structure
        The slab that should be on the bottom.
    d : float, optional
        The desired distance between the slabs. The default is 2.5.

    Returns
    -------
    t_copy : pymatgen.core.structure.Structure
        Top slab that is shifted so that the lowest atom is at +d/2
    b_copy : pymatgen.core.structure.Structure
        Bottom slab that is shifted so that the topmost atom is at -d/2

    """
    t_copy = top_slab.copy()
    b_copy = bottom_slab.copy()
    top_zs = []
    bot_zs = []
    for s in t_copy.sites:
        top_zs.append(s.coords[-1])
    top_shift = -min(top_zs) + d / 2

    for s in b_copy.sites:
        bot_zs.append(s.coords[-1])
    bot_shift = -max(bot_zs) - d / 2

    t_copy.translate_sites(
        indices=range(len(t_copy.sites)),
        vector=[0, 0, top_shift],
        frac_coords=False,
        to_unit_cell=False,
    )
    b_copy.translate_sites(
        indices=range(len(b_copy.sites)),
        vector=[0, 0, bot_shift],
        frac_coords=False,
        to_unit_cell=False,
    )
    return t_copy, b_copy


def interface_name(
        mpid1: str,
        mpid2: str,
        miller1: tuple,
        miller2: tuple,
        shift1: float = None,
        shift2: float = None,
) -> str:
    """Return a name for an interface based on MP-IDs and miller indices.

    Parameters
    ----------
    mpid1 : str
        MP-ID of the first material.
    mpid2 : str
        MP-ID of the second material.
    miller1 : list of int, or str
        Miller indices of the first material given either as list or str with
        3 letters.
    miller2 : list of int, or str
        Miller indices of the second material given either as list or str with
        3 letters.
    shift1 : float, optional
        Shift of the first slab that determines its terminations. The default is None.
    shift2 : float, optional
        Shift of the second slab that determines its terminations. The default is None.

    Returns
    -------
    name : str
        Unique name for the interface of two slabs.

    """

    mp_connection = MPConnection()
    f1 = mp_connection.get_property_from_mp(
        mpid=mpid1, properties=["formula_pretty"]
    )
    f1 = f1["formula_pretty"]

    f2 = mp_connection.get_property_from_mp(
        mpid=mpid2, properties=["formula_pretty"]
    )
    f2 = f2["formula_pretty"]

    if type(miller1) is list:
        m1 = "".join(str(s) for s in miller1)
    else:
        m1 = miller1
    if type(miller2) is list:
        m2 = "".join(str(s) for s in miller2)
    else:
        m2 = miller2

    if shift1 is None:
        name1 = f"{f1 + m1}-({mpid1})"
    else:
        name1 = f"{f1 + m1}_{shift1}-({mpid1})"

    if shift2 is None:
        name2 = f"{f2 + m2}-({mpid2})"
    else:
        name2 = f"{f2 + m2}_{shift2}-({mpid2})"

    name = "_".join(sorted([name1, name2]))
    return name


def make_pymatgen_slab(
        bulk_struct: Structure,
        miller: list,
        min_thickness: int = 8,
        min_vacuum: float = 20,
) -> Slab:
    """
    Uses pmg slab generator to create 1 slab from a bulk struct and miller index.

    First a conventional unit cell is created using SpacegroupAnalyzer, then
    possible magnetic moments get transferred and the slab is created using
    pymatgen.core.surface.SlabGenerator. Only the first slab of the returned
    list is taken! Oxidation states are guessed, so polarity of the slab might
    be queried.

    Parameters
    ----------
    bulk_struct : pymatgen.core.structure.Structure
        the bulk structure from which the slab will be constructed.
    miller : list
        A list of three miller indices
    min_thickness : int, optional
        Thickness of the slab in layers (probably will end up thicker).
        The default is 8.
    min_vacuum : float, optional
        Thickness of the vacuum layer in Angstrom. The default is 20.

    Returns
    -------
    pymatgen.core.surface.Slab
        pymatgen Slab object.

    """

    bulk_conv = SpacegroupAnalyzer(
        bulk_struct
    ).get_conventional_standard_structure(keep_site_properties=True)
    bulk_conv = transfer_average_magmoms(bulk_struct, bulk_conv)

    sg = SlabGenerator(
        initial_structure=bulk_conv,
        miller_index=miller,
        center_slab=True,
        primitive=True,
        lll_reduce=True,
        max_normal_search=max([abs(m) for m in miller]),
        min_slab_size=min_thickness,
        min_vacuum_size=min_vacuum,
    )

    slab = sg.get_slabs(
        bonds=None,
        ftol=0.1,
        tol=0.1,
        max_broken_bonds=0,
        symmetrize=False,
        repair=False,
    )[0]
    slab.add_oxidation_state_by_guess()
    return slab


def flip_slab(slab: Slab) -> Slab:
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
    flip_matrix = np.array(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]
    )
    flipped_coords = np.dot(slab.cart_coords, flip_matrix)

    try:
        flipped_slab = Slab(
            lattice=slab.lattice,
            species=slab.species,
            coords=flipped_coords,
            miller_index=slab.miller_index,
            oriented_unit_cell=slab.oriented_unit_cell,
            shift=slab.shift,
            scale_factor=slab.scale_factor,
            reconstruction=slab.reconstruction,
            coords_are_cartesian=True,
            site_properties=slab.site_properties,
        )
    except:
        flipped_slab = Structure(
            lattice=slab.lattice,
            species=slab.species,
            coords=flipped_coords,
            coords_are_cartesian=True,
            site_properties=slab.site_properties,
        )
    return center_slab(flipped_slab)
