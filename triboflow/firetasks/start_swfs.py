from atomate.utils.utils import env_chk
from fireworks import explicit_serialize, FiretaskBase, FWAction
from pymatgen.core import Structure
from pymatgen.core.interface import Interface
from pymatgen.core.surface import Slab

from hitmen_utils.db_tools import VaspDB
from triboflow.workflows.subworkflows import (
    adhesion_energy_swf,
    calc_pes_swf,
    calc_ppes_swf,
    converge_swf,
    charge_analysis_swf,
)


@explicit_serialize
class FT_StartChargeAnalysisSWF(FiretaskBase):
    """Start a charge redistribution analysis subworkflow.

    Take an interface from the high_level database and compute
    the charge density redistribution through a subworkflow.

    Parameters
    ----------
    mp_id_1 : str
        MaterialsProject ID number for the first material
    mp_id_2 : str
        MaterialsProject ID number for the second material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    external_pressure : float
        External pressure in GPa.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    high_level : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    interface_label : str, optional
        Label that the relaxed interface has in the high-level database. Can
        be either structure@min (default), or structure@max at the moment.
    """

    required_params = [
        "mp_id_1",
        "mp_id_2",
        "functional",
        "external_pressure",
    ]
    optional_params = ["db_file", "high_level", "interface_label"]

    def run_task(self, fw_spec):
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)
        interface_label = self.get("interface_label", "relaxed_structure@min")

        db = VaspDB(db_file, high_level=hl_db)

        name = fw_spec.get("interface_name")

        interface_dict = db.find_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
        )

        redistribution_was_calculated = interface_dict.get(
            "charge_density_redist"
        )
        comp_params = interface_dict.get("comp_parameters", {})

        if not redistribution_was_calculated:
            interface = Interface.from_dict(interface_dict[interface_label])

            SWF = charge_analysis_swf(
                interface=interface,
                interface_name=name,
                functional=functional,
                external_pressure=external_pressure,
                db_file=db_file,
                high_level=hl_db,
                comp_parameters=comp_params,
            )

            return FWAction(detours=SWF)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartAdhesionSWF(FiretaskBase):
    """Start an adhesion subworkflow.

    Take relaxed top and bottom slabs of an interface, as well as the relaxed
    interface structure (by default the one with the lowest energy) and compute
    the adhesion energy through a subworkflow.

    Parameters
    ----------
    mp_id_1 : str
        MaterialsProject ID number for the first material
    mp_id_2 : str
        MaterialsProject ID number for the second material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    external_pressure : float
        External pressure in GPa.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    adhesion_handle: str, optional
        Flag under which the adhesion energy will be saved in the
        interface_data collection of the high_level database.
    high_level : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    """

    required_params = [
        "mp_id_1",
        "mp_id_2",
        "functional",
        "external_pressure",
    ]
    optional_params = ["db_file", "adhesion_handle", "high_level"]

    def run_task(self, fw_spec):
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        adhesion_handle = self.get("adhesion_handle", "adhesion_energy@min")
        hl_db = self.get("high_level", True)

        db = VaspDB(db_file, high_level=hl_db)

        name = fw_spec.get("interface_name")

        interface_dict = db.find_data(
            collection=functional + ".interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
        )

        adhesion_was_calculated = interface_dict.get(adhesion_handle)
        comp_params = interface_dict.get("comp_parameters", {})

        if not adhesion_was_calculated:
            top_slab = Slab.from_dict(interface_dict["top_aligned_relaxed"])
            bottom_slab = Slab.from_dict(
                interface_dict["bottom_aligned_relaxed"]
            )
            interface = Structure.from_dict(
                interface_dict["relaxed_structure@min"]
            )

            swf = adhesion_energy_swf(
                top_slab,
                bottom_slab,
                interface,
                external_pressure=external_pressure,
                interface_name=name,
                functional=functional,
                comp_parameters=comp_params,
                db_file=db_file,
                high_level=hl_db,
            )

            return FWAction(detours=swf)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartBulkConvoSWF(FiretaskBase):
    """Starts a convergence subworkflow.

    Starts either an energy cutoff or kpoint density convergence of a material
    with a given MPID and functional through a subworkflow.

    Parameters
    ----------
    conv_type : str
        Either "kpoints" or "encut", depending on what to converge.
    mp_id : str
        MaterialsProject ID number for the material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    encut_start : float, optional
        Starting encut value for the first run. Defaults to the largest EMIN
        in the POTCAR.
    encut_incr : float, optional
        Increment for the encut during the convergence. Defaults to 25.
    k_dens_start : float, optional
        Starting kpoint density in 1/Angstrom. Defaults to 1.0
    k_dens_incr : float, optional
        Increment for the kpoint convergence. Can be set quite small since
        there is a check in place to see if a new mesh is actually constructed
        for each density. Defaults to 0.1.
    n_converge : int, optional
        Number of calculations that have to be inside the convergence
        threshold for convergence to be reached. Defaults to 3.
    high_level : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    """

    _fw_name = "Start Encut or Kdensity Convergence"
    required_params = ["conv_type", "mp_id", "functional"]
    optional_params = [
        "db_file",
        "encut_start",
        "encut_incr",
        "k_dens_start",
        "k_dens_incr",
        "n_converge",
        "high_level",
    ]

    def run_task(self, fw_spec):
        conv_type = self.get("conv_type")
        mp_id = self.get("mp_id")
        functional = self.get("functional")
        db_file = self.get("db_file")
        n_converge = self.get("n_converge", 3)
        encut_start = self.get("encut_start", None)
        encut_incr = self.get("encut_incr", 25)
        k_dens_start = self.get("k_dens_start", 2.0)
        k_dens_incr = self.get("k_dens_incr", 0.1)
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)

        if conv_type not in ["kpoints", "encut"]:
            raise ValueError(
                '"type" input must be either "kpoints" or'
                '"encut".\nYou have passed {}'.format(conv_type)
            )

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        data = db_high.find_data(
            collection=f"{functional}.bulk_data", fltr={"mpid": mp_id}
        )

        if conv_type == "encut":
            stop_convergence = data.get("encut_info")
        elif conv_type == "kpoints":
            stop_convergence = data.get("k_dense_info")

        if not stop_convergence:
            structure_dict = data.get("structure_equiVol")
            if not structure_dict:
                structure_dict = data.get("primitive_structure")
                if not structure_dict:
                    structure_dict = data.get("structure_fromMP")
                    if not structure_dict:
                        raise LookupError(
                            "No structure found that can be used "
                            "as input for the convergence swf."
                        )
            structure = Structure.from_dict(structure_dict)
            comp_params = data.get("comp_parameters", {})
            SWF = converge_swf(
                structure=structure,
                conv_type=conv_type,
                flag=mp_id,
                comp_parameters=comp_params,
                functional=functional,
                encut_start=encut_start,
                encut_incr=encut_incr,
                k_dens_start=k_dens_start,
                k_dens_incr=k_dens_incr,
                n_converge=n_converge,
                print_help=False,
            )
            return FWAction(detours=SWF, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartPESCalcSWF(FiretaskBase):
    """Start a PES subworkflow.

    Starts a PES subworkflow using data from the high-level database.
    This is intended to be used to start a PES subworkflow from a main
    workflow.

    Parameters
    ----------
    mp_id_1 : str
        Materials Project database ID for the first material of the interface.
    mp_id_2 : str
        Materials Project database ID for the second material of the interface.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    external_pressure : float
        External pressure in GPa.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    prerelax : bool, optional
        Whether to perform a prerelaxation using a network potential
        before starting a DFT relaxation. Defaults to True.
    prerelax_calculator : str, optional
        Which network potential to use for the prerelaxation. Defaults
        to 'm3gnet'.
    prerelax_kwargs : dict, optional
        Keyword arguments to be passed to the ASE calculator for
        the prerelaxation.

    Returns
    -------
    FWAction that produces a detour PES subworkflow.
    """

    required_params = [
        "mp_id_1",
        "mp_id_2",
        "functional",
        "external_pressure",
    ]
    optional_params = [
        "db_file",
        "high_level",
        "prerelax",
        "prerelax_calculator",
        "prerelax_kwargs",
    ]

    def run_task(self, fw_spec):
        functional = self.get("functional")
        external_pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)
        prerelax = self.get("prerelax", True)
        prerelax_calculator = self.get("prerelax_calculator", "m3gnet")
        prerelax_kwargs = self.get("prerelax_kwargs", {})

        name = fw_spec.get("interface_name")

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        interface_dict = db_high.find_data(
            collection=f"{functional}.interface_data",
            fltr={"name": name, "external_pressure": external_pressure},
        )
        comp_params = interface_dict["comp_parameters"]
        interface = Interface.from_dict(interface_dict["unrelaxed_structure"])
        already_done = interface_dict.get("relaxed_structure@min")

        if not already_done:
            SWF = calc_pes_swf(
                interface=interface,
                interface_name=name,
                functional=functional,
                external_pressure=external_pressure,
                comp_parameters=comp_params,
                output_dir=None,
                prerelax=prerelax,
                prerelax_calculator=prerelax_calculator,
                prerelax_kwargs=prerelax_kwargs,
                db_file=db_file,
                high_level=hl_db,
            )

            return FWAction(detours=SWF, update_spec=fw_spec)

        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartPPESWF(FiretaskBase):
    """
    Start a CalcPPES_SWF subworkflow that calculates a PPES.

    The workflow is only added if there are not already relevant results in
    the high-level database.

    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    distance_list : list of float, optional
        Modification of the equilibrium distance between the slabs.
        The default is [-0.5, -0.25, 0.0, 0.25, 0.5, 2.5, 3.0, 4.0, 5.0, 7.5].
    out_name : str, optional
        Name for the PPES data in the high-level database. The default is
        'PPES@minimum'.
    structure_name : str, optional
        Name of the structure in the interface entry to the high-level database
        for which the PPES should be calculated. The default is
        'minimum_relaxed'.
    spec : dict, optional
        fw_spec that can be passed to the SWF and will be passed on. The
        default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        Subworkflow to calculate the PPES for a certain interface.

    """

    required_params = ["interface_name", "functional", "distance_list"]
    optional_params = [
        "db_file",
        "structure_name",
        "out_name",
        "high_level",
    ]

    def run_task(self, fw_spec):
        name = self.get("interface_name")
        functional = self.get("functional")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)

        structure_name = self.get("structure_name", "minimum_relaxed")
        out_name = self.get("out_name", "PPES@minimum")

        d_list = self.get("distance_list")
        hl_db = self.get("high_level", True)

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        interface_dict = db_high.find_data(
            collection=f"{functional}.interface_data",
            fltr={"name": name, "external_pressure": 0.0},
        )
        if not interface_dict:
            raise ValueError(
                f"Interface {name} with functional {functional} and "
                f"0 external_pressure not found in high-level database {db_file}!"
            )

        calc_PPES = True
        if interface_dict.get("PPES") is not None:
            if interface_dict["PPES"].get(out_name) is not None:
                print(
                    "\n A PPES-object with out_name: "
                    + out_name
                    + "\n has already been created in the interface entry: "
                    + name
                    + "\n for the "
                    + functional
                    + " functional."
                )
                calc_PPES = False

        if calc_PPES:
            SWF = calc_ppes_swf(
                interface_name=name,
                functional=functional,
                distance_list=d_list,
                out_name=out_name,
                structure_name=structure_name,
                spec=fw_spec,
            )

            return FWAction(additions=SWF, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)
