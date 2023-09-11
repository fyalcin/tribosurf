#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python


from fireworks import Firework

from triboflow.firetasks.PES import (
    FT_FindHighSymmPoints,
    FT_StartPESCalcs,
    FT_RetrievePESEnergies,
    FT_ComputePES,
)
from triboflow.utils.file_manipulation import copy_output_files

__author__ = "Michael Wolloch"
__copyright__ = "Copyright 2020, Michael Wolloch"
__version__ = "0.1"
__maintainer__ = "Michael Wolloch"
__email__ = "michael.wolloch@univie.ac.at"
__date__ = "March 11th, 2020"


# =============================================================================
# Custom FireWorks
# =============================================================================


def run_pes_calc_fw(
    interface,
    interface_name,
    external_pressure,
    functional,
    comp_parameters,
    tag,
    FW_name,
    prerelax=True,
    prerelax_calculator="m3gnet",
    prerelax_kwargs=None,
    db_file="auto",
    high_level=True,
):
    """Compute high-symmetry points for an interface and start PES calculations.

    Combines two Fireworks that find the high-symmetry points for the interface
    and start the VASP calculations for the unique high-symmetry points
    respectively.

    :param interface: Interface object for which the PES is to be constructed
    :type interface: pymatgen.core.interface.Interface

    :param interface_name: Unique name for the interface that is used in the output and the database.
    :type interface_name: str

    :param external_pressure: External pressure in GPa.
    :type external_pressure: float

    :param functional: Functional to be used. 'PBE' or 'SCAN' will work.
    :type functional: str

    :param comp_parameters: Dictionary containing computational options. E.g. encut, k_dens, vdw,...
    :type comp_parameters: dict

    :param tag: combination of the interface_name and a uuid. To uniquely identify the computations in the database.
    :type tag: str

    :param FW_name: Name of the Firework.
    :type FW_name: str

    :param prerelax: Whether to perform a prerelaxation using a neural-network potential before starting a DFT relaxation. Defaults to True.
    :type prerelax: bool, optional

    :param prerelax_calculator: Which neural-network potential to use for the prerelaxation. Defaults to 'm3gnet'.
    :type prerelax_calculator: str, optional

    :param prerelax_kwargs: Keyword arguments to be passed to the ASE calculator for the prerelaxation.
    :type prerelax_kwargs: dict, optional

    :param db_file: Full path to the db.json file that should be used. Defaults to '>>db_file<<', to use env_chk.
    :type db_file: str

    :param high_level: Whether to use the high-level database. Defaults to True.
    :type high_level: bool, optional

    :return: First Firework of a PES subworkflow.
    :rtype: fireworks.core.firework.Firework

    """
    if prerelax_kwargs is None:
        prerelax_kwargs = {}
    ft_1 = FT_FindHighSymmPoints(
        interface=interface,
        interface_name=interface_name,
        functional=functional,
        external_pressure=external_pressure,
        db_file=db_file,
        high_level=high_level,
    )

    ft_2 = FT_StartPESCalcs(
        interface_name=interface_name,
        comp_parameters=comp_parameters,
        tag=tag,
        external_pressure=external_pressure,
        prerelax=prerelax,
        prerelax_calculator=prerelax_calculator,
        prerelax_kwargs=prerelax_kwargs,
        db_file=db_file,
    )

    fw = Firework([ft_1, ft_2], name=FW_name)

    return fw


def make_pes_fw(
    interface_name,
    functional,
    external_pressure,
    tag,
    FW_name,
    file_output,
    output_dir,
    remote_copy=False,
    server=None,
    user=None,
    port=None,
    db_file="auto",
    high_level=True,
):
    """Retrieve PES calculations from the database and compute the PES.

    Retrieve the computed energies of the unique high-symmetry points and match
    them to the replicate points. Duplicates the points, interpolates with
    radial basis functions and saves the results. Plots the results as well.
    Optionally write file output and copy it to an output directory.

    :param interface_name: Unique name for the interface that is used in the output and the database.
    :type interface_name: str

    :param functional: Functional to be used. 'PBE' or 'SCAN' will work.
    :type functional: str

    :param external_pressure: External pressure in GPa.
    :type external_pressure: float

    :param tag: combination of the interface_name and a uuid. To uniquely identify the computations in the database.
    :type tag: str

    :param FW_name: Name of the Firework.
    :type FW_name: str

    :param file_output: Determines if files are written to disk.
    :type file_output: bool

    :param output_dir: Location the output files are copied to if file_output is selected.
    :type output_dir: str

    :param remote_copy: If true, scp will be used to copy the results to a remote server. Be advised that ssh-key certification must be set up between the two machines. The default is False.
    :type remote_copy: bool, optional

    :param server: Fully qualified domain name of the server the output should be copied to. The default is None.
    :type server: str, optional

    :param user: The username on the remote server.
    :type user: str, optional

    :param port: On some machines ssh-key certification is only supported for certain ports. A port may be selected here. The default is None.
    :type port: int, optional

    :param db_file: Full path to the db.json file that should be used. Defaults to '>>db_file<<', to use env_chk.
    :type db_file: str

    :param high_level: Whether to use the high-level database. Defaults to True.
    :type high_level: bool, optional

    :return: Final Firework of a PES subworkflow.
    :rtype: fireworks.core.firework.Firework

    """
    ft_1 = FT_RetrievePESEnergies(
        interface_name=interface_name,
        functional=functional,
        tag=tag,
        external_pressure=external_pressure,
        db_file=db_file,
        high_level=high_level,
    )
    ft_2 = FT_ComputePES(
        interface_name=interface_name,
        functional=functional,
        external_pressure=external_pressure,
        file_output=file_output,
        db_file=db_file,
        high_level=high_level,
    )

    if file_output:
        output_files = [str(interface_name) + ".png"]
        ft_3 = copy_output_files(
            file_list=output_files,
            output_dir=output_dir,
            remote_copy=remote_copy,
            server=server,
            user=user,
            port=port,
        )

        fw = Firework([ft_1, ft_2, ft_3], name=FW_name)
    else:
        fw = Firework([ft_1, ft_2], name=FW_name)

    return fw
