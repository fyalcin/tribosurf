"""Firetasks for converging the energy cutof in the triboflow project.

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
from fireworks import FWAction, FiretaskBase, Firework, Workflow, FileWriteTask
from fireworks.utilities.fw_utilities import explicit_serialize

from hitmen_utils.kpoints import MeshFromDensity
from hitmen_utils.vasp_tools import (
    get_custom_vasp_static_settings,
    get_emin_and_emax,
)
from hitmen_utils.db_tools import VaspDB
from triboflow.utils.file_manipulation import copy_output_files
from triboflow.utils.utils import is_list_converged


@explicit_serialize
class FT_UpdateBMLists(FiretaskBase):
    """Fetch information about the EOS fit from the DB and update the lists.

    Used with FT_Convo to converge energy cutoffs or kpoints using the
    get_wf_bulk_modulus workflow of atomate. This Firetasks reads the
    neccessary information from the eos collection of the database. Since no
    useful tag is placed, the identification of the correct entry is done
    by the chemical formula and the timestamp. The shared data entry for the
    convergence is then identified via a tag and updated with the new
    equilibrium volume and the bulk modulus.

    Parameters
    ----------
    formula : str
        Chemical formula on the material to be matched with the database.
    tag : str
        String from an uuid4 to identify the shared data entry in the database.
    db_file : str, optional
        Full path to the db.json file detailing access to the database.
        Defaults to '>>db_file<<' to use with env_chk.

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
        results = db.vasp_calc_db.db.eos.find(
            {"formula_pretty": formula}
        ).sort("created_at", pymongo.DESCENDING)
        BM = results["bulk_modulus"]
        V0 = results["results"]["v0"]

        # Update data arrays in the database
        db = VaspDB(db_file=db_file)
        db.update_data(
            collection="BM_data_sharing",
            fltr={"tag": tag},
            new_values={"$push": {"BM_list": BM, "V0_list": V0}},
        )


@explicit_serialize
class FT_Convo(FiretaskBase):
    """Converge the encut or kpoint density via fits to an EOS.

    Uses the get_bulk_modulus workflow of atomate to fit Birch-Murnaghen EOS
    for increasing values of the energy cutoff or kpoint density. Once bulk
    modulus and equilibrium volume are converged, the subsequent detours are
    stopped and the convergence data is passed on.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-pint grids.
    conv_type : str
        Either "kpoints" or "encut", depending on what to converge.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
    tag : str
        String from an uuid4 to identify the shared data entry in the database
        and group the deformations calculations together.
    deformations: list of lists, optional
        List of deformation matrices for the fit to the EOS. Defaults to None,
        which results in 5 volumes from 90% to 110% of the initial volume.
    n_converge : int, optional
        Number of iterations that have to be below the convergence threshold
        for the system to be considered converged. Defaults to 3
    encut_start : float, optional
        Starting encut value for the first run in eV. Defaults to the largest
        EMIN in the POTCAR.
    encut_incr : float, optional
        Increment for the encut during the convergence. Defaults to 25.
    k_dens_start : float, optional
        Starting kpoint density in 1/Angstrom. Defaults to 1.0
    k_dens_increment : float, optional
        Increment for the kpoint convergence. Can be set quite small since
        there is a check in place to see if a new mesh is actually constructed
        for each density. Defaults to 0.1.
    k_dens_default : float, optional
        Default (quite high) kpoints density for encut convergence studies if
        no k_dens parameter is found in the comp_parameters. The default is 12.5
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The username on the remote server. The default is None.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.

    Returns
    -------
    FWActions that produce detour subworkflows until convergence is reached.
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

        V0_tolerance = comp_params.get("volume_tolerence", 0.001)
        BM_tolerance = comp_params.get("BM_tolerence", 0.01)

        # Get the data arrays from the database (returns None when not there)
        db = VaspDB(db_file=db_file)
        data = nav.find_data(collection="BM_data_sharing", fltr={"tag": tag})

        if data:
            BM_list = data.get("BM_list")
            V0_list = data.get("V0_list")
            convo_list = data.get("convo_list")
        else:
            BM_list = None
            V0_list = None
            convo_list = None

        if BM_list is None:
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

            BM_WF = get_wf_bulk_modulus(
                struct,
                deformations,
                vasp_input_set=vis,
                vasp_cmd=VASP_CMD,
                db_file=db_file,
                eos="birch_murnaghan",
                tag=tag,
            )

            formula = struct.composition.reduced_formula
            UAL_FW = Firework(
                [
                    FT_UpdateBMLists(formula=formula, tag=tag),
                    FT_Convo(
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

            BM_WF.append_wf(Workflow.from_Firework(UAL_FW), BM_WF.leaf_fw_ids)
            # Use add_modify_incar powerup to add KPAR and NCORE settings
            # based on env_chk in my_fworker.yaml
            BM_WF = add_modify_incar(BM_WF)
            # Set up the entry for the data arrays in the database
            set_data = {
                "tag": tag,
                "chem_formula": formula,
                "created_on": str(datetime.now()),
                "convo_list": convo_list,
                "BM_list": [],
                "V0_list": [],
            }
            nav.insert_data(collection="BM_data_sharing", data=set_data)

            return FWAction(detours=BM_WF)

        else:
            BM_tol = BM_list[-1] * BM_tolerance
            V0_tol = V0_list[-1] * V0_tolerance

            if is_list_converged(
                BM_list, BM_tol, n_converge
            ) and is_list_converged(V0_list, V0_tol, n_converge):
                # Handle the last iteration
                final_BM = BM_list[-n_converge]
                final_V0 = V0_list[-n_converge]
                flag = self.get("flag")
                functional = self.get("functional")

                scaled_structure = struct.copy()
                scaled_structure.scale_lattice(final_V0)
                struct_dict = scaled_structure.as_dict()

                if conv_type == "encut":
                    final_encut = convo_list[-n_converge]
                    print("")
                    print(" Convergence reached for BM and cell volume.")
                    print(
                        " Final encut = {} eV; Final BM = {} GPa;"
                        "Final Volume = {} Angstrom³".format(
                            final_encut, final_BM, final_V0
                        )
                    )
                    print("")
                    print(" The scaled output structure is:\n")
                    pprint(struct_dict)
                    output_dict = {
                        "encut_info": {
                            "final_encut": final_encut,
                            "final_BM": final_BM,
                            "final_volume": final_V0,
                            "BM_list": BM_list,
                            "V0_list": V0_list,
                            "convo_list": convo_list,
                            "BM_tol_abs": BM_tol,
                            "BM_tol_rel": BM_tolerance,
                            "V0_tol_abs": V0_tol,
                            "V0_tol_rel": V0_tolerance,
                        },
                        "equilibrium_volume": final_V0,
                        "bulk_moduls": final_BM,
                        "comp_parameters.encut": final_encut,
                        "structure_equiVol": struct_dict,
                    }
                else:
                    final_k_dens = convo_list[-n_converge]
                    print("")
                    print(" Convergence reached for total energy per atom.")
                    print(
                        " Final k_dens = {}; Final BM = {} GPa;"
                        "Final Volume = {} Angstrom³".format(
                            final_k_dens, final_BM, final_V0
                        )
                    )
                    print("")
                    print("")
                    output_dict = {
                        "k_dense_info": {
                            "final_k_dens": final_k_dens,
                            "final_BM": final_BM,
                            "final_volume": final_V0,
                            "BM_list": BM_list,
                            "V0_list": V0_list,
                            "convo_list": convo_list,
                            "BM_tol_abs": BM_tol,
                            "BM_tol_rel": BM_tolerance,
                            "V0_tol_abs": V0_tol,
                            "V0_tol_rel": V0_tolerance,
                        },
                        "equilibrium_volume": final_V0,
                        "bulk_moduls": final_BM,
                        "comp_parameters.k_dens": final_k_dens,
                        "structure_equiVol": struct_dict,
                    }

                db_high = VaspDB(db_file=db_file, high_level=hl_db)
                db_high.update_data(
                    collection=functional + ".bulk_data",
                    fltr={"mpid": flag},
                    new_values={"$set": output_dict},
                    upsert=True,
                )

                if conv_type == "encut":
                    nav.update_data(
                        collection="BM_data_sharing",
                        fltr={"tag": tag},
                        new_values={
                            "$set": {
                                "final_encut": final_encut,
                                "final_BM": final_BM,
                                "final_volume": final_V0,
                                "BM_tol_abs": BM_tol,
                                "BM_tol_rel": BM_tolerance,
                                "V0_tol_abs": V0_tol,
                                "V0_tol_rel": V0_tolerance,
                            }
                        },
                    )
                else:
                    nav.update_data(
                        collection="BM_data_sharing",
                        fltr={"tag": tag},
                        new_values={
                            "$set": {
                                "final_k_dense": final_k_dens,
                                "final_BM": final_BM,
                                "final_volume": final_V0,
                                "BM_tol_abs": BM_tol,
                                "BM_tol_rel": BM_tolerance,
                                "V0_tol_abs": V0_tol,
                                "V0_tol_rel": V0_tolerance,
                            }
                        },
                    )

                # handle file output:
                if file_output:
                    write_FT = FileWriteTask(
                        files_to_write=[
                            {
                                "filename": flag + "_output_dict.txt",
                                "contents": pformat(output_dict),
                            }
                        ]
                    )

                    copy_FT = copy_output_files(
                        file_list=[flag + "_output_dict.txt"],
                        output_dir=output_dir,
                        remote_copy=remote_copy,
                        server=server,
                        user=server,
                        port=port,
                    )

                    FW = Firework(
                        [write_FT, copy_FT],
                        name="Copy Convergence SWF results",
                    )

                    WF = Workflow.from_Firework(
                        FW, name="Copy Convergence SWF results"
                    )

                    return FWAction(update_spec=fw_spec, detours=WF)
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
                KPTS = MeshFromDensity(
                    struct, k_dens, compare_density=convo_list[-1]
                )
                while KPTS.are_meshes_the_same():
                    k_dens = k_dens + k_dens_incr
                    KPTS = MeshFromDensity(
                        struct, k_dens, compare_density=convo_list[-1]
                    )
                # Pass kspacing to ensure correct meshes for all deformations
                comp_params["kspacing"] = 1.0 / k_dens

            vis = get_custom_vasp_static_settings(
                struct, comp_params, "bulk_from_scratch"
            )

            BM_WF = get_wf_bulk_modulus(
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
            UAL_FW = Firework(
                [
                    FT_UpdateBMLists(formula=formula, tag=tag),
                    FT_Convo(
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

            BM_WF.append_wf(Workflow.from_Firework(UAL_FW), BM_WF.leaf_fw_ids)
            # Use add_modify_incar powerup to add KPAR and NCORE settings
            # based on env_chk in my_fworker.yaml
            BM_WF = add_modify_incar(BM_WF)

            # Update Database entry for convo list
            if conv_type == "encut":
                nav.update_data(
                    collection="BM_data_sharing",
                    fltr={"tag": tag},
                    new_values={"$push": {"convo_list": encut}},
                )
            else:
                nav.update_data(
                    collection="BM_data_sharing",
                    fltr={"tag": tag},
                    new_values={"$push": {"convo_list": k_dens}},
                )

            return FWAction(detours=BM_WF)
