"""Firetasks for converging the energy cutoff in the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""
from datetime import datetime
from pprint import pprint, pformat

import numpy as np
import pymongo
from atomate.utils.utils import env_chk
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.workflows.base.bulk_modulus import get_wf_bulk_modulus
from fireworks import (
    FWAction,
    FiretaskBase,
    Firework,
    Workflow,
    FileWriteTask,
)
from fireworks.utilities.fw_utilities import explicit_serialize

from htflow_utils.db_tools import VaspDB
from htflow_utils.kpoints import MeshFromDensity
from htflow_utils.vasp_tools import (
    get_custom_vasp_static_settings,
    get_emin_and_emax,
)
from triboflow.utils.file_manipulation import copy_output_files


@explicit_serialize
class UpdateBMLists(FiretaskBase):
    """Fetch information about the EOS fit from the DB and update the lists.

    Used with Converge to converge energy cutoffs or kpoints using the
    get_wf_bulk_modulus workflow of atomate. This Firetasks reads the
    necessary information from the eos collection of the database. Since no
    useful tag is placed, the identification of the correct entry is done
    by the chemical formula and the timestamp. The shared data entry for the
    convergence is then identified via a tag and updated with the new
    equilibrium volume and the bulk modulus.

    :param formula: Chemical formula on the material to be matched with the
        database.
    :type formula: str

    :param tag: String from an uuid4 to identify the shared data entry in the
        database.
    :type tag: str

    :param db_file: Full path to the db.json file detailing access to the
        database. Defaults to '>>db_file<<' to use with env_chk.
    :type db_file: str, optional

    :return: None
    :rtype: NoneType
    """

    _fw_name = "Update Bulk Modulus Lists"
    required_params = ["formula", "tag"]
    optional_params = ["db_file"]

    def run_task(self, fw_spec):
        formula = self.get("formula")
        tag = self.get("tag")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)

        db = VaspDB(db_file=db_file)
        results = list(
            db.vasp_calc_db.db.eos.find({"formula_pretty": formula}).sort(
                "created_at", pymongo.DESCENDING
            )
        )[0]
        bm = results["bulk_modulus"]
        v0 = results["results"]["v0"]

        # Update data arrays in the database
        db = VaspDB(db_file=db_file)
        db.update_data(
            collection="BM_data_sharing",
            fltr={"tag": tag},
            new_values={"$push": {"BM_list": bm, "V0_list": v0}},
        )


@explicit_serialize
class Converge(FiretaskBase):
    """Converge the encut or kpoint density via fits to an EOS.

    Uses the get_bulk_modulus workflow of atomate to fit Birch-Murnaghan EOS
    for increasing values of the energy cutoff or kpoint density. Once bulk
    modulus and equilibrium volume are converged, the subsequent detours are
    stopped and the convergence data is passed on.

    :param structure: The structure for which to converge the K-point grids or the
        energy cutoff.
    :type structure: pymatgen.core.structure.Structure

    :param conv_type: Either "kpoints" or "encut", depending on what to
        converge.
    :type conv_type: str

    :param comp_params: Dictionary of computational parameters for the VASP
        calculations.
    :type comp_params: dict

    :param tag: String from an uuid4 to identify the shared data entry in the
        database and group the deformations calculations together.
    :type tag: str

    :param deformations: List of deformation matrices for the fit to the EOS.
        Defaults to None, which results in 5 volumes from 90% to 110% of the
        initial volume.
    :type deformations: list of lists, optional

    :param n_converge: Number of iterations that have to be below the
        convergence threshold for the system to be considered converged.
        Defaults to 3
    :type n_converge: int, optional

    :param encut_start: Starting encut value for the first run in eV.
        Defaults to the largest EMIN in the POTCAR.
    :type encut_start: float, optional

    :param encut_incr: Increment for the encut during the convergence.
        Defaults to 25.
    :type encut_incr: float, optional

    :param k_dens_start: Starting kpoint density in 1/Angstrom.
        Defaults to 1.0
    :type k_dens_start: float, optional

    :param k_dens_increment: Increment for the kpoint convergence.
        Can be set quite small since there is a check in place to see if a new
        mesh is actually constructed for each density. Defaults to 0.1.
    :type k_dens_increment: float, optional

    :param k_dens_default: Default (quite high) kpoints density for encut
        convergence studies if no k_dens parameter is found in the
        comp_parameters. The default is 12.5
    :type k_dens_default: float, optional

    :param db_file: Full path to the db.json file that should be used.
        Defaults to '>>db_file<<', to use env_chk.
    :type db_file: str

    :param file_output: Toggles file output. The default is False.
    :type file_output: bool, optional

    :param output_dir: Defines a directory the output is to be copied to.
        (Do not use a trailing / and/or relative location symbols like ~/.)
        The default is None.
    :type output_dir: str, optional

    :param remote_copy: If true, scp will be used to copy the results to a
        remote server. Be advised that ssh-key certification must be set up
        between the two machines. The default is False.
    :type remote_copy: bool, optional

    :param server: Fully qualified domain name of the server the output should
        be copied to. The default is None.
    :type server: str, optional

    :param user: The username on the remote server. The default is None.
    :type user: str, optional

    :param port: On some machines ssh-key certification is only supported for
        certain ports. A port may be selected here. The default is None.
    :type port: int, optional

    :return: FWActions that produce detour subworkflows until convergence is
        reached.
    :rtype: FWAction
    """

    _fw_name = "Convergence for encut or kpoints"
    required_params = [
        "structure",
        "conv_type",
        "comp_params",
        "tag",
        "flag",
        "functional",
    ]
    optional_params = [
        "deformations",
        "n_converge",
        "encut_start",
        "encut_incr",
        "k_dens_start",
        "k_dens_incr",
        "k_dens_default",
        "db_file",
        "file_output",
        "output_dir",
        "remote_copy",
        "server",
        "user",
        "port",
        "high_level",
    ]

    def run_task(self, fw_spec):
        deforms = []
        for i in np.arange(0.95, 1.05, 0.025):
            dm = np.eye(3) * i
            deforms.append(dm)
        n_converge = self.get("n_converge", 3)
        encut_start = self.get("encut_start", None)
        encut_incr = self.get("encut_incr", 25)
        k_dens_start = self.get("k_dens_start", 2.0)
        k_dens_incr = self.get("k_dens_incr", 0.1)
        k_dens_def = self.get("k_dens_default", 12.5)
        deformations = self.get("deformations")
        file_output = self.get("file_output", False)
        output_dir = self.get("output_dir", None)
        remote_copy = self.get("remote_copy", False)
        server = self.get("server", None)
        user = self.get("user", None)
        port = self.get("port", None)
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)
        if not deformations:
            deformations = deforms

        struct = self["structure"]
        conv_type = self.get("conv_type")
        comp_params = self["comp_params"]
        tag = self["tag"]

        v0_tolerance = comp_params.get("volume_tolerance", 0.001)
        bm_tolerance = comp_params.get("BM_tolerance", 0.01)

        # Get the data arrays from the database (returns None when not there)
        db = VaspDB(db_file=db_file)
        data = db.find_data(collection="BM_data_sharing", fltr={"tag": tag})

        if data:
            bm_list = data.get("BM_list")
            v0_list = data.get("V0_list")
            convo_list = data.get("convo_list")
        else:
            bm_list = None
            v0_list = None
            convo_list = None

        if bm_list is None:
            # Handle the first iteration.

            if conv_type == "encut":
                if not encut_start:
                    vis = get_custom_vasp_static_settings(
                        struct, comp_params, "bulk_from_scratch"
                    )
                    # Get the largest EMIN value of the potcar and round up to the
                    # next whole 25.
                    encut_dict = get_emin_and_emax(vis.potcar)
                    enmax = encut_dict["ENMAX"]
                    encut_start = int(25 * np.ceil(enmax / 25))
                comp_params["encut"] = encut_start
                convo_list = [encut_start]
                # Pass kspacing to ensure correct meshes for all deformations
                if "k_dens" in comp_params:
                    comp_params["kspacing"] = 1.0 / comp_params["k_dens"]
                else:
                    comp_params["kspacing"] = 1.0 / k_dens_def
            else:
                convo_list = [k_dens_start]
                # Pass kspacing to ensure correct meshes for all deformations
                comp_params["kspacing"] = 1.0 / k_dens_start

            vis = get_custom_vasp_static_settings(
                struct, comp_params, "bulk_from_scratch"
            )

            bm_wf = get_wf_bulk_modulus(
                struct,
                deformations,
                vasp_input_set=vis,
                vasp_cmd=VASP_CMD,
                db_file=db_file,
                eos="birch_murnaghan",
                tag=tag,
            )

            formula = struct.composition.reduced_formula
            ual_fw = Firework(
                [
                    UpdateBMLists(formula=formula, tag=tag),
                    Converge(
                        structure=struct,
                        conv_type=conv_type,
                        comp_params=comp_params,
                        tag=tag,
                        flag=self["flag"],
                        functional=self["functional"],
                        db_file=db_file,
                        encut_incr=encut_incr,
                        encut_start=encut_start,
                        k_dens_start=k_dens_start,
                        k_dens_incr=k_dens_incr,
                        file_output=file_output,
                        output_dir=output_dir,
                        remote_copy=remote_copy,
                        server=server,
                        user=user,
                        port=port,
                        high_level=hl_db,
                    ),
                ],
                name="Update BM Lists and Loop",
            )

            bm_wf.append_wf(Workflow.from_Firework(ual_fw), bm_wf.leaf_fw_ids)
            # Use add_modify_incar powerup to add KPAR and NCORE settings
            # based on env_chk in my_fworker.yaml
            bm_wf = add_modify_incar(bm_wf)
            # Set up the entry for the data arrays in the database
            set_data = {
                "tag": tag,
                "chem_formula": formula,
                "created_on": str(datetime.now()),
                "convo_list": convo_list,
                "BM_list": [],
                "V0_list": [],
            }
            db.insert_data(collection="BM_data_sharing", data=set_data)

            return FWAction(detours=bm_wf)

        else:
            bm_tol = bm_list[-1] * bm_tolerance
            v0_tol = v0_list[-1] * v0_tolerance

            if is_list_converged(
                bm_list, bm_tol, n_converge
            ) and is_list_converged(v0_list, v0_tol, n_converge):
                # Handle the last iteration
                final_bm = bm_list[-n_converge]
                final_v0 = v0_list[-n_converge]
                flag = self.get("flag")
                functional = self.get("functional")

                scaled_structure = struct.copy()
                scaled_structure.scale_lattice(final_v0)
                struct_dict = scaled_structure.as_dict()

                if conv_type == "encut":
                    final_encut = convo_list[-n_converge]
                    print("")
                    print(" Convergence reached for BM and cell volume.")
                    print(
                        " Final encut = {} eV; Final BM = {} GPa;"
                        "Final Volume = {} Angstrom³".format(
                            final_encut, final_bm, final_v0
                        )
                    )
                    print("")
                    print(" The scaled output structure is:\n")
                    pprint(struct_dict)
                    output_dict = {
                        "encut_info": {
                            "final_encut": final_encut,
                            "final_BM": final_bm,
                            "final_volume": final_v0,
                            "BM_list": bm_list,
                            "V0_list": v0_list,
                            "convo_list": convo_list,
                            "BM_tol_abs": bm_tol,
                            "BM_tol_rel": bm_tolerance,
                            "V0_tol_abs": v0_tol,
                            "V0_tol_rel": v0_tolerance,
                        },
                        "equilibrium_volume": final_v0,
                        "bulk_moduls": final_bm,
                        "comp_parameters.encut": final_encut,
                        "structure_equiVol": struct_dict,
                    }
                elif conv_type == "kpoints":
                    final_k_dens = convo_list[-n_converge]
                    print("")
                    print(" Convergence reached for total energy per atom.")
                    print(
                        " Final k_dens = {}; Final BM = {} GPa;"
                        "Final Volume = {} Angstrom³".format(
                            final_k_dens, final_bm, final_v0
                        )
                    )
                    print("")
                    print("")
                    output_dict = {
                        "k_dense_info": {
                            "final_k_dens": final_k_dens,
                            "final_BM": final_bm,
                            "final_volume": final_v0,
                            "BM_list": bm_list,
                            "V0_list": v0_list,
                            "convo_list": convo_list,
                            "BM_tol_abs": bm_tol,
                            "BM_tol_rel": bm_tolerance,
                            "V0_tol_abs": v0_tol,
                            "V0_tol_rel": v0_tolerance,
                        },
                        "equilibrium_volume": final_v0,
                        "bulk_moduls": final_bm,
                        "comp_parameters.k_dens": final_k_dens,
                        "structure_equiVol": struct_dict,
                    }
                else:
                    raise ValueError(
                        "Convergence type {} not recognized!".format(conv_type)
                    )

                db_high = VaspDB(db_file=db_file, high_level=hl_db)
                db_high.update_data(
                    collection=functional + ".bulk_data",
                    fltr={"mpid": flag},
                    new_values={"$set": output_dict},
                    upsert=True,
                )

                if conv_type == "encut":
                    db.update_data(
                        collection="BM_data_sharing",
                        fltr={"tag": tag},
                        new_values={
                            "$set": {
                                "final_encut": final_encut,
                                "final_BM": final_bm,
                                "final_volume": final_v0,
                                "BM_tol_abs": bm_tol,
                                "BM_tol_rel": bm_tolerance,
                                "V0_tol_abs": v0_tol,
                                "V0_tol_rel": v0_tolerance,
                            }
                        },
                    )
                else:
                    db.update_data(
                        collection="BM_data_sharing",
                        fltr={"tag": tag},
                        new_values={
                            "$set": {
                                "final_k_dense": final_k_dens,
                                "final_BM": final_bm,
                                "final_volume": final_v0,
                                "BM_tol_abs": bm_tol,
                                "BM_tol_rel": bm_tolerance,
                                "V0_tol_abs": v0_tol,
                                "V0_tol_rel": v0_tolerance,
                            }
                        },
                    )

                # handle file output:
                if file_output:
                    write_ft = FileWriteTask(
                        files_to_write=[
                            {
                                "filename": flag + "_output_dict.txt",
                                "contents": pformat(output_dict),
                            }
                        ]
                    )

                    copy_ft = copy_output_files(
                        file_list=[flag + "_output_dict.txt"],
                        output_dir=output_dir,
                        remote_copy=remote_copy,
                        server=server,
                        user=server,
                        port=port,
                    )

                    fw = Firework(
                        [write_ft, copy_ft],
                        name="Copy Convergence SWF results",
                    )

                    wf = Workflow.from_Firework(
                        fw, name="Copy Convergence SWF results"
                    )

                    return FWAction(update_spec=fw_spec, detours=wf)
                else:
                    return FWAction(update_spec=fw_spec)

            # Make a normal iteration

            if conv_type == "encut":
                encut = convo_list[-1] + encut_incr
                comp_params["encut"] = encut
                # Pass kspacing to ensure correct meshes for all deformations
                if "k_dens" in comp_params:
                    comp_params["kspacing"] = 1.0 / comp_params["k_dens"]
                else:
                    comp_params["kspacing"] = 1.0 / k_dens_def
            else:
                k_dens = convo_list[-1] + k_dens_incr
                # Ensure that the new density leads to a different mesh.
                kpts = MeshFromDensity(
                    struct, k_dens, compare_density=convo_list[-1]
                )
                while kpts.are_meshes_the_same():
                    k_dens = k_dens + k_dens_incr
                    kpts = MeshFromDensity(
                        struct, k_dens, compare_density=convo_list[-1]
                    )
                # Pass kspacing to ensure correct meshes for all deformations
                comp_params["kspacing"] = 1.0 / k_dens

            vis = get_custom_vasp_static_settings(
                struct, comp_params, "bulk_from_scratch"
            )

            bm_wf = get_wf_bulk_modulus(
                struct,
                deformations,
                vasp_input_set=vis,
                vasp_cmd=VASP_CMD,
                db_file=DB_FILE,
                # user_kpoints_settings=uks,
                eos="birch_murnaghan",
                tag=tag,
            )

            formula = struct.composition.reduced_formula
            ual_fw = Firework(
                [
                    UpdateBMLists(formula=formula, tag=tag),
                    Converge(
                        structure=struct,
                        conv_type=conv_type,
                        comp_params=comp_params,
                        tag=tag,
                        flag=self["flag"],
                        functional=self["functional"],
                        db_file=db_file,
                        encut_incr=encut_incr,
                        encut_start=encut_start,
                        k_dens_start=k_dens_start,
                        k_dens_incr=k_dens_incr,
                        file_output=file_output,
                        output_dir=output_dir,
                        remote_copy=remote_copy,
                        server=server,
                        user=user,
                        port=port,
                    ),
                ],
                name="Update BM Lists and Loop",
            )

            bm_wf.append_wf(Workflow.from_Firework(ual_fw), bm_wf.leaf_fw_ids)
            # Use add_modify_incar powerup to add KPAR and NCORE settings
            # based on env_chk in my_fworker.yaml
            bm_wf = add_modify_incar(bm_wf)

            # Update Database entry for convo list
            if conv_type == "encut":
                db.update_data(
                    collection="BM_data_sharing",
                    fltr={"tag": tag},
                    new_values={"$push": {"convo_list": encut}},
                )
            else:
                db.update_data(
                    collection="BM_data_sharing",
                    fltr={"tag": tag},
                    new_values={"$push": {"convo_list": k_dens}},
                )

            return FWAction(detours=bm_wf)


def is_list_converged(input_list: list[float], tol: float, n: int = 3) -> bool:
    """Check if the last n values of an array are within tol of each other.

    :param input_list: List of values to be checked for convergence.
    :type input_list: list of float

    :param tol: Tolerance for the convergence.
    :type tol: float

    :param n: Number of entries at the end of input_list that have to be within
        tol for the list to be considered converged. The default is 3.
    :type n: int, optional

    :return: True if input_list is converged, False otherwise.
    :rtype: bool
    """
    if len(input_list) <= n:
        return False
    else:
        check_list = [False] * n
        tmp = input_list.copy()
        tmp.reverse()
        for i, b in enumerate(check_list):
            if abs(tmp[0] - tmp[i + 1]) < tol:
                check_list[i] = True
        return all(check_list)
