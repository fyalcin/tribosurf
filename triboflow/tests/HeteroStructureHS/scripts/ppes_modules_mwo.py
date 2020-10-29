from pymatgen.core.sites import PeriodicSite


def StackAlignedSlabs(bottom_slab, top_slab, top_shift=[0,0,0]):
    """
    Combine slabs that are centered around 0 into a single structure.
    
    Optionally shift the top slab by a vector of cartesian coordinates.

    Parameters
    ----------
    bottom_slab : pymatgen.core.structure.Structure
        Bottom slab.
    top_slab : pymatgen.core.structure.Structure
        Top slab.
    top_shift : list of 3 floats, optional
        Vector of caresian coordinates with which to shift the top slab.
        The default is [0,0,0].

    Returns
    -------
    interface : pymatgen.core.structure.Structure
        DESCRIPTION.

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
                                species=s.species)
        interface.sites.append(new_site)
    
    return interface

def ReCenterAlignedSlabs(top_slab, bottom_slab, d=2.5):
    """
    Center two slabs around z=0 and give them the distance d.

    Parameters
    ----------
    top_slab : pymatgen.core.structboture.Structure
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

