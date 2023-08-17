import subprocess

from pymatgen.io.vasp.inputs import Kpoints

from triboflow.utils.file_manipulation import remove_matching_files


def get_generalized_kmesh(structure, k_dist, RemoveSymm=False, Vasp6=True):
    """Get a generalized Monkhorst Pack mesh for a given structure.

    Prepares the necessary files (POSCAR, PRECALC, and, if the structure has
    a 'magmom' site property also INCAR) for the K-Point Grid Generator of the
    Mueller group at John Hopkins http://muellergroup.jhu.edu/K-Points.html
    Runs the getKPoints script and reads the KPOINT file produced into a
    pymatgen Kpoints object.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        Structure for which the kpoints grid is to be generated.
    k_dist : float
        The minimum allowed distance between lattice points on the real-space
        superlattice. This determines the density of the k-point grid. A larger
        value will result in denser grids.

    Returns
    -------
    KPTS : pymatgen.io.vasp.inputs.Kpoints
        Pymatgen Kpoints object representing a generalized MP mesh.

    """

    precalc = [
        "MINDISTANCE = " + str(k_dist),
        "MINTOTALKPOINTS = 4",
        "GAPDISTANCE = 6 ",
        "MONOCLINIC_SEARCH_DEPTH = 2500",
        "TRICLINIC_SEARCH_DEPTH = 1500",
    ]
    if Vasp6:
        precalc.append("WRITE_LATTICE_VECTORS = True")

    if RemoveSymm in ["STRUCTURAL", "TIME_REVERSAL", "ALL"]:
        precalc.append("REMOVE_SYMMETRY = " + RemoveSymm)
    with open("PRECALC", "w") as out:
        for line in precalc:
            out.write(line + "\n")

    magmom_list = structure.site_properties.get("magmom")
    if magmom_list:
        with open("INCAR", "w") as out:
            out.write("ISPIN = 2")
            out.write("MAGMOM = " + " ".join(str(m) for m in magmom_list))

    structure.to(fmt="poscar", filename="POSCAR")
    get_kpoints_file = subprocess.Popen("getKPoints")
    get_kpoints_file.communicate()
    KPTS = Kpoints().from_file("KPOINTS")
    remove_matching_files(["KPOINTS*", "POSCAR*", "INCAR", "PRECALC"])

    return KPTS
