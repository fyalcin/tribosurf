import warnings

import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, Slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.ext.matproj import MPRester
from pymatgen.transformations.standard_transformations import (
    DeformStructureTransformation,
)

from hitmen_utils.shaper import Shaper
from triboflow.utils.database import  Navigator
from triboflow.utils.mp_connection import MPConnection



def get_interface_distance(structure):
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


def get_SG_from_mpid(
    mpid,
    functional,
    miller,
    db_file="auto",
    high_level=True,
    min_slab=15,
    min_vac=15,
    lll_reduce=True,
    center=True,
    in_unit_planes=True,
    prim=True,
    reorient_lattice=True,
):
    """
    Generates the SlabGenerator for the given material, miller index and other parameters.

    Parameters
    ----------
    mpid : str
        Unique MaterialsProject ID describing the structure.
    functional : str
        Functional used for the bulk calculation. Used for accessing
        the correct collection in the database to load the conventional
        standard structure.
    miller : list
        Desired orientation of the SlabGenerator object.
    db_file : str, optional
        Full path of the db.json file. If 'auto', looks it up from the
        FW_CONFIG_FILE environment variable. The default is 'auto'.
    high_level : str, bool optional
        High level database in which to search for the primitive bulk
        structure. If True, looks it up from the db.json file.
        The default is True.
    min_slab : float or int, optional
        Minimum thickness of the slab region in either Angstroms or number
        of layers depending on in_unit_planes. The default is 15.
    min_vac : float or int, optional
        Minimum thickness of the slab region in either Angstroms or number
        of layers depending on in_unit_planes. The default is 15.
    lll_reduce : bool, optional
        Whether to perform LLL orthogonalization or not on the slabs generated
        by the SlabGenerator.
        The default is True.
    center : bool, optional
        Whether to center the slabs generated by the SlabGenerator or not.
        The default is True.
    in_unit_planes : bool, optional
        Decides whether min_slab and min_vac describe the thicknesses of the
        corresponding regions in number of layers or Angstroms. If 'True',
        number of planes is used. The default is True.
    prim : bool, optional
        Whether to look for smaller oriented unit cells after initializing
        the SlabGenerator. Leads to smaller slabs. The default is True.
    reorient_lattice : bool, optional
        Whether to orient the c vector of the slabs in the original z direction.
        The default is True.

    Returns
    -------
    SG : pymatgen.core.surface.SlabGenerator
        Pymatgen SlabGenerator object with the given parameters.

    """
    coll = f"{functional}.bulk_data"
    bulk_conv = get_conv_bulk_from_mpid(mpid, coll, db_file, high_level)
    max_normal_search = max([abs(m) for m in miller])
    SG = SlabGenerator(
        bulk_conv,
        miller,
        min_slab,
        min_vac,
        lll_reduce,
        center,
        in_unit_planes,
        prim,
        max_normal_search,
        reorient_lattice,
    )
    return SG


def get_conv_bulk_from_mpid(mpid, coll, db_file="auto", high_level=True):
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
    nav = Navigator(db_file, high_level)
    bulk_dict = nav.find_data(coll, {"mpid": mpid})
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
    ).get_conventional_standard_structure()
    return bulk_conv


def slab_from_file(
    filename, mpid, functional, miller, db_file="auto", high_level=True
):
    """
    Loads up a surface from a file (can be CIF, VASP, POSCAR ..) and turn it
    into a pymatgen Slab object.

    Parameters
    ----------
    filename : str
        Full path of the file to be loaded.
    mpid : str
        Unique MaterialsProject ID describing the structure.
    functional : str
        Functional used in the bulk calculation. Used to generate
        the correct SlabGenerator.
    miller : list
        Desired orientation of the SlabGenerator object to be used in the
        generation of the Slab.
    db_file : str, optional
        Full path of the db.json file. If 'auto', looks it up from the
        FW_CONFIG_FILE environment variable. The default is 'auto'.
    high_level : str, bool optional
        High level database in which to search for the primitive bulk
        structure. If True, looks it up from the db.json file.
        The default is True.

    Returns
    -------
    slab : pymatgen.core.surface.Slab
        Pymatgen Slab object having the same parameters as the structure given
        in the file. shift, scale_factor, and in some cases the oriented_unit_cell
        will not match the structure in the file.

    """
    struct = Structure.from_file(filename)
    SG = get_SG_from_mpid(mpid, functional, miller, db_file, high_level)
    tmp_slab = SG.get_slabs()[0]
    ouc = Shaper.get_matching_ouc(tmp_slab)
    if ouc:
        scale_factor = struct.lattice.a / ouc.lattice.a
        scale_array = [
            (scale_factor, 0, 0),
            (0, scale_factor, 0),
            (0, 0, scale_factor),
        ]
        deformer = DeformStructureTransformation(deformation=scale_array)
        ouc = deformer.apply_transformation(ouc)
    else:
        print(
            f"Matching OUC could not be found for the structure loaded from file. A non-matching\n"
            f"OUC will be used instead. Results are still useful."
        )
        ouc = SG.oriented_unit_cell
    slab = Slab(
        lattice=struct.lattice,
        species=struct.species,
        coords=struct.frac_coords,
        miller_index=miller,
        oriented_unit_cell=ouc,
        shift=0,
        scale_factor=SG.slab_scale_factor,
    )
    return slab, SG


def transfer_average_magmoms(magnetic_struct, struct_without_magmoms):
    """Set magmom for a structure based on the average value of each species of a reference structure.

    For unit cells of the same structure, it is not always trivial to transfer
    the site properties. This function attempts to transfer at least the magmom
    site property between two structures with the same species, but not
    necessarily the same number of sites. For each species the average value
    of the magentic moments in the magnetic input structure is computed and
    set as a site property for all atoms of the same species in the output
    structure. NOTE THAT THIS WILL GIVE GENERALLY WRONG RESULTS FOR ALL BUT
    SIMPLE FERROMAGENTIC STRUCTURES!

    Parameters
    ----------
    magnetic_struct : pymatgen.core.structure.Structure
        Input structure with "magmom" site property.
    struct_without_magmoms : pymatgen.core.structure.Structure
        Input structure with no "magmom" site property but the same species.

    Returns
    -------
    new_struct : pymatgen.core.structure.Structure
        copy of struct_without_magmoms with added "magmom" site property.

    """

    mag_struct = magnetic_struct.copy()
    new_struct = struct_without_magmoms.copy()

    if not mag_struct.site_properties.get("magmom"):
        print("No magnetic moments to transfer. Doing nothing...")
        return new_struct

    if not sorted(mag_struct.types_of_species) == sorted(
        new_struct.types_of_species
    ):
        warnings.warn(
            "\n##################################################\n"
            "You are trying to transfer magnetig moments between\n"
            "two structures which contain different species and\n"
            "                 THIS CANNOT WORK!\n"
            "The code will continue to run, without transfering\n"
            "any magnetic moments. Convergence might be slow..."
            "\n##################################################\n"
        )
        return new_struct

    magmom_dict = {}
    for s in mag_struct.types_of_species:
        magmom_dict[s] = []
        for i, el in enumerate(mag_struct.species):
            if s == el:
                magmom_dict[s].append(
                    mag_struct.site_properties.get("magmom")[i]
                )
        magmom_dict[s] = np.mean(magmom_dict[s])

    new_magmoms = []
    for s in new_struct.species:
        new_magmoms.append(magmom_dict[s])
    new_struct.add_site_property("magmom", new_magmoms)

    return new_struct


def slab_from_structure(miller, structure):
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


def clean_up_site_properties(structure):
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
                if m == None:
                    new_magmom.append(0.0)
                else:
                    new_magmom.append(m)
            struct.add_site_property("magmom", new_magmom)
        else:
            if any(struct.site_properties[key]) == None:
                struct.remove_site_property(key)
    return struct


def stack_aligned_slabs(bottom_slab, top_slab, top_shift=[0, 0, 0]):
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
        Vector of caresian coordinates with which to shift the top slab.
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


def recenter_aligned_slabs(top_slab, bottom_slab, d=2.5):
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


def interface_name(mp_id_1, miller_1, mp_id_2, miller_2):
    """Return a name for an interface based on MP-IDs and miller indices.

    Parameters
    ----------
    mp_id_1 : str
        MP-ID of the first material.
    miller_1 : list of int, or str
        Miller indices of the first material given either as list or str with
        3 letters.
    mp_id_2 : str
        MP-ID of the second material.
    miller_2 : list of int, or str
        Miller indices of the second material given either as list or str with
        3 letters.

    Returns
    -------
    name : str
        Unique name for the interface of two slabs.

    """

    mp_connection = MPConnection()
    f1 = mp_connection.get_property_from_mp(
        mp_id=mp_id_1, properties=["formula_pretty"]
    )
    f1 = f1["formula_pretty"]

    f2 = mp_connection.get_property_from_mp(
        mp_id=mp_id_2, properties=["formula_pretty"]
    )
    f2 = f2["formula_pretty"]

    if type(miller_1) is list:
        m1 = "".join(str(s) for s in miller_1)
    else:
        m1 = miller_1
    if type(miller_2) is list:
        m2 = "".join(str(s) for s in miller_2)
    else:
        m2 = miller_2

    n1 = min(f1 + m1, f2 + m2)
    n2 = max(f1 + m1, f2 + m2)
    ids = min(mp_id_1 + "_" + mp_id_2, mp_id_2 + "_" + mp_id_1)
    name = "_".join((n1, n2, ids))
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
    possible magnetic moments get transfered and the slab is created using
    pymatgen.core.surface.SlabGenerator. Only the first slab of the returned
    list is taken! Oxidation states are guessed, so polarity of the slab might
    be quaried.

    Parameters
    ----------
    bulk_struct : pymatgen.core.structure.Structure
        the bulk structure from which the slab will be constructed.
    miller : list
        List of three miller indices
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
    ).get_conventional_standard_structure()
    bulk_conv = transfer_average_magmoms(bulk_struct, bulk_conv)

    SG = SlabGenerator(
        initial_structure=bulk_conv,
        miller_index=miller,
        center_slab=True,
        primitive=True,
        lll_reduce=True,
        max_normal_search=max([abs(l) for l in miller]),
        min_slab_size=min_thickness,
        min_vacuum_size=min_vacuum,
    )

    slab = SG.get_slabs(
        bonds=None,
        ftol=0.1,
        tol=0.1,
        max_broken_bonds=0,
        symmetrize=False,
        repair=False,
    )[0]
    slab.add_oxidation_state_by_guess()
    return slab
