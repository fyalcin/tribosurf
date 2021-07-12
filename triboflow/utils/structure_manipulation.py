import numpy as np
import warnings
from pymatgen.core.surface import Slab
from pymatgen.core.sites import PeriodicSite

from triboflow.utils.database import NavigatorMP


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
    
    if not mag_struct.site_properties.get('magmom'):
        print('No magnetic moments to transfer. Doing nothing...')
        return new_struct
    
    if not sorted(mag_struct.types_of_species) == sorted(new_struct.types_of_species):
        warnings.warn('\n##################################################\n'
                      'You are trying to transfer magnetig moments between\n'
                      'two structures which contain different species and\n'
                      '                 THIS CANNOT WORK!\n'
                      'The code will continue to run, without transfering\n'
                      'any magnetic moments. Convergence might be slow...'
                      '\n##################################################\n')
        return new_struct
    
    magmom_dict={}
    for s in mag_struct.types_of_species:
        magmom_dict[s]=[]
        for i, el in enumerate(mag_struct.species):
            if s==el:
                magmom_dict[s].append(mag_struct.site_properties.get('magmom')[i])
        magmom_dict[s] = np.mean(magmom_dict[s])
        
    new_magmoms = []
    for s in new_struct.species:
        new_magmoms.append(magmom_dict[s])
    new_struct.add_site_property('magmom', new_magmoms)
    
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
    return Slab(lattice=structure.lattice,
                species=structure.species_and_occu,
                coords=structure.frac_coords,
                miller_index=miller,
                oriented_unit_cell=structure,
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=structure.site_properties)

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
        if key == 'magmom':
            new_magmom = []
            for m in struct.site_properties[key]:
                if m == None:
                    new_magmom.append(0.0)
                else:
                    new_magmom.append(m)
            struct.add_site_property('magmom', new_magmom)
        else:
            struct.remove_site_property(key)
    return struct

def stack_aligned_slabs(bottom_slab, top_slab, top_shift=[0,0,0]):
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
    
    t_copy.translate_sites(indices=range(len(t_copy.sites)),
                           vector=top_shift,
                           frac_coords=False, to_unit_cell=False)
    
    for s in t_copy.sites:
        new_site = PeriodicSite(lattice=interface.lattice,
                                coords=s.frac_coords,
                                coords_are_cartesian=False,
                                species=s.species,
                                properties=s.properties)
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
    top_zs=[]
    bot_zs=[]
    for s in t_copy.sites:
        top_zs.append(s.coords[-1])
    top_shift = -min(top_zs) + d/2
    
    for s in b_copy.sites:
        bot_zs.append(s.coords[-1])
    bot_shift = -max(bot_zs) - d/2

    t_copy.translate_sites(indices=range(len(t_copy.sites)),
                           vector=[0,0,top_shift],
                           frac_coords=False, to_unit_cell=False)
    b_copy.translate_sites(indices=range(len(b_copy.sites)),
                           vector=[0,0,bot_shift],
                           frac_coords=False, to_unit_cell=False)
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

    nav_mp = NavigatorMP()
    f1 = nav_mp.get_property_from_mp(
        mp_id=mp_id_1,
        properties=['pretty_formula'])
    f1 = f1['pretty_formula']

    f2 = nav_mp.get_property_from_mp(
        mp_id=mp_id_2,
        properties=['pretty_formula'])
    f2 = f2['pretty_formula']

    if type(miller_1) is list:
        m1 = ''.join(str(s) for s in miller_1)
    else:
        m1 = miller_1
    if type(miller_2) is list:
        m2 = ''.join(str(s) for s in miller_2)
    else:
        m2 = miller_2

    n1 = min(f1+m1, f2+m2)
    n2 = max(f1+m1, f2+m2)
    ids = min(mp_id_1+'_'+mp_id_2, mp_id_2+'_'+mp_id_1)
    name = '_'.join((n1, n2, ids))
    return name
