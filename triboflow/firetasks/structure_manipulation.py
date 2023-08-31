""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""
import itertools
import numpy as np
from atomate.utils.utils import env_chk
from fireworks import FWAction, FiretaskBase, Firework, Workflow, FileWriteTask
from fireworks.utilities.fw_utilities import explicit_serialize
from pprint import pprint, pformat
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from uuid import uuid4

from hitmen_utils.db_tools import VaspDB
from hitmen_utils.vasp_tools import get_custom_vasp_relax_settings
from hitmen_utils.workflows import dynamic_relax_swf
from surfen.utils.structure_manipulation import add_bulk_to_db
from triboflow.phys.interface_matcher import (
    InterfaceMatcher,
    get_consolidated_comp_params,
)
from triboflow.utils.file_manipulation import copy_output_files
from triboflow.utils.structure_manipulation import (
    interface_name,
    transfer_average_magmoms,
)


@explicit_serialize
class FT_StartBulkPreRelax(FiretaskBase):
    """Start a subworkflow as a detour to relax the cell shape and positions
    of a primitive structure depending on lattice parameters, and then move
    the optimized primitive structure to the high level database.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    encut : float, optional
        Energy cutoff for the relaxation run. Defaults to 1000.
    k_dens : int, optional
        kpoint density in 1/Angstrom. Defaults to a (quite high) 15.
    high_level : str, optional
        Name of the high level database the structure should be queried in
        and later the results written to. Defaults to 'triboflow'.
    Returns
    -------
        Starts a subworkflow as a detour to the current workflow.
    """

    _fw_name = "Start a cell shape relaxation"
    required_params = ["mp_id", "functional"]
    optional_params = ["db_file", "encut", "k_dens", "high_level"]

    def run_task(self, fw_spec):
        mp_id = self.get("mp_id")
        functional = self.get("functional")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)

        # Querying the structure from the high level database.
        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        data = db_high.find_data(
            collection=f"{functional}.bulk_data", fltr={"mpid": mp_id}
        )

        prim_struct = data.get("primitive_structure")
        if not prim_struct:
            struct = data.get("structure_fromMP")
            if not struct:
                raise LookupError(
                    "No structure found in the database that can "
                    "be used as input for cell shape relaxation."
                )
            struct = Structure.from_dict(struct)
            prim_struct = SpacegroupAnalyzer(struct).get_primitive_standard_structure()
            prim_struct = transfer_average_magmoms(struct, prim_struct)
        else:
            prim_struct = Structure.from_dict(prim_struct)

        a = np.round(prim_struct.lattice.a, 6)
        b = np.round(prim_struct.lattice.b, 6)
        c = np.round(prim_struct.lattice.c, 6)

        if data.get("pre_relaxed") or ((a == c) and (b == c)):
            return FWAction(update_spec=fw_spec)
        else:
            # Querying the computational parameters from the high level
            # database and updating with the optional inputs
            comp_params = data.get("comp_parameters")
            encut = self.get("encut", 1000)
            k_dens = self.get("k_dens", 15)
            comp_params.update({"encut": encut, "k_dens": k_dens})

            tag = "CellShapeRelax-{}".format(str(uuid4()))
            vis = get_custom_vasp_relax_settings(
                prim_struct, comp_params, "bulk_pos_shape_relax"
            )
            RelaxWF = dynamic_relax_swf(
                inputs_list=[[prim_struct, vis, tag]],
                prerelax_system=True,
                prerelax_kwargs={"relax_cell": True},
            )

            MoveResultsFW = Firework(
                [
                    FT_UpdatePrimStruct(
                        functional=functional,
                        tag=tag,
                        flag=mp_id,
                        high_level=hl_db,
                    )
                ],
                name="Move pre-relaxed structure for {}".format(prim_struct.formula),
            )
            MoveResultsWF = Workflow([MoveResultsFW])

            RelaxWF.append_wf(MoveResultsWF, RelaxWF.leaf_fw_ids)

            return FWAction(detours=RelaxWF, update_spec=fw_spec)


@explicit_serialize
class FT_UpdatePrimStruct(FiretaskBase):
    """Update the primitive structure in the high level database.

    Parameters
    ----------
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    flag : str
        An identifier to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise, use
        something unique which you can find again.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    Returns
    -------
        Moves the optimized primitive structure from the low level
        database to the high level database.
    """

    _fw_name = "Update primitive structure in the high level DB"
    required_params = ["functional", "tag", "flag"]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        functional = self.get("functional")
        tag = self.get("tag")
        flag = self.get("flag")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level", True)

        db = VaspDB(db_file=db_file)
        calc = db.find_data("tasks", {"task_label": tag})
        out = calc["output"]

        struct_dict = {
            "primitive_structure": out["structure"],
            "pre_relaxed": True,
        }

        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        db_high.update_data(
            collection=functional + ".bulk_data",
            fltr={"mpid": flag},
            new_values={"$set": struct_dict},
            upsert=False,
        )

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_GetRelaxedSlab(FiretaskBase):
    """Get the relaxed structure, and put a Slab into the high-level DB.

    Parameters
    ----------
    flag : str
        An identifier to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise, use
        something unique which you can find again.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    tag : str
        Unique identifier for the Optimize FW.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    struct_out_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'relaxed_slab'.
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
        The username on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.

    Returns
    -------
        Pymatgen Slab into the high-level DB.
    """

    required_params = ["flag", "miller", "functional", "tag"]
    optional_params = [
        "db_file",
        "struct_out_name",
        "file_output",
        "high_level",
        "output_dir",
        "remote_copy",
        "server",
        "user",
        "port",
    ]

    def run_task(self, fw_spec):
        flag = self.get("flag")
        if type(self["miller"]) == str:
            miller = [int(k) for k in list(self["miller"])]
        else:
            miller = self["miller"]

        functional = self.get("functional")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)

        hl_db = self.get("high_level", True)
        out_name = self.get("struct_out_name", "relaxed_slab")
        file_output = self.get("file_output", False)
        output_dir = self.get("output_dir", None)
        remote_copy = self.get("remote_copy", False)
        server = self.get("server", None)
        user = self.get("user", None)
        port = self.get("port", None)

        # Check if a relaxed slab is already in the DB entry
        db_high = VaspDB(db_file=db_file, high_level=hl_db)
        slab_data = db_high.find_data(
            collection=f"{functional}.slab_data",
            fltr={"mpid": flag, "miller": miller},
        )

        if out_name not in slab_data:
            # Get results from OptimizeFW
            db = VaspDB(db_file=db_file)
            vasp_calc = db.find_data(
                collection="tasks", fltr={"task_label": self["tag"]}
            )
            relaxed_slab = Structure.from_dict(vasp_calc["output"]["structure"])
            slab = Slab(
                relaxed_slab.lattice,
                relaxed_slab.species_and_occu,
                relaxed_slab.frac_coords,
                miller,
                Structure.from_sites(relaxed_slab, to_unit_cell=True),
                shift=0,
                scale_factor=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                site_properties=relaxed_slab.site_properties,
            )

            db_high = VaspDB(db_file=db_file, high_level=hl_db)
            db_high.update_data(
                collection=functional + ".slab_data",
                fltr={"mpid": flag, "miller": miller},
                new_values={"$set": {out_name: slab.as_dict()}},
            )
        else:
            db = VaspDB(db_file=db_file)
            vasp_calc = db.find_data(
                collection="tasks", fltr={"task_label": self["tag"]}
            )

            if vasp_calc:
                relaxed_slab = Structure.from_dict(vasp_calc["output"]["structure"])
                slab = Slab(
                    relaxed_slab.lattice,
                    relaxed_slab.species_and_occu,
                    relaxed_slab.frac_coords,
                    miller,
                    Structure.from_sites(relaxed_slab, to_unit_cell=True),
                    shift=0,
                    scale_factor=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    site_properties=relaxed_slab.site_properties,
                )
                print("")
                print(
                    " A slab with the selected output name already exists in the DB."
                    " It will not be overwritten with the new relaxed slab.\n"
                    " If needed you can update the data manually."
                )
                print("")
            else:
                slab = Slab.from_dict(slab_data[out_name])
                print("")
                print(
                    " A slab with the selected output name already exists in the DB."
                    " No new slab has been relaxed.\n"
                )

        # screen output:
        print("")
        print("Relaxed output structure as pymatgen.surface.Slab dictionary:")
        pprint(slab.as_dict())
        print("")

        # handle file output:
        if file_output:
            poscar_str = Poscar(slab).get_string()
            poscar_name = flag + "_Relaxed_slab_POSCAR.vasp"
            slab_name = flag + "_Relaxed_slab_dict.txt"
            write_FT = FileWriteTask(
                files_to_write=[
                    {"filename": poscar_name, "contents": poscar_str},
                    {
                        "filename": slab_name,
                        "contents": pformat(slab.as_dict()),
                    },
                ]
            )
            copy_FT = copy_output_files(
                file_list=[poscar_name, slab_name],
                output_dir=output_dir,
                remote_copy=remote_copy,
                server=server,
                user=user,
                port=port,
            )
            FW = Firework([write_FT, copy_FT], name="Copy SlabRelax SWF results")
            WF = Workflow.from_Firework(FW, name="Copy SlabRelax SWF results")

            return FWAction(update_spec=fw_spec, detours=WF)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_AddBulkToDB(FiretaskBase):
    _fw_name = "Add bulk to DB"
    required_params = ["mpid", "functional"]
    optional_params = ["db_file", "high_level", "custom_data"]

    def run_task(self, fw_spec):
        mpid = self.get("mpid")
        functional = self.get("functional")
        db_file = self.get("db_file", "auto")
        high_level = self.get("high_level", True)
        custom_data = self.get("custom_data", None)

        add_bulk_to_db(
            mpid=mpid,
            coll=f"{functional}.bulk_data",
            db_file=db_file,
            high_level=high_level,
            custom_data=custom_data,
        )

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_MakeHeteroStructure(FiretaskBase):
    """Matches two slab systems to form a heterogeneous interface.

    If the match fails, the accuracy criteria are relaxed in steps of 5% until
    a match is found.

    Parameters
    ----------
    mp_id_1 : str
        ID number from the MP to identify material 1.
    mp_id_2 : str
        ID number from the MP for material 2.
    functional : str
        functional that is used for the calculation.
    external_pressure : float
        External pressure in GPa.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.

    Returns
    -------
    Matched interface structure and bottom and top slabs in interface_data DB.
    """

    _fw_name = "Make Hetero Structure"
    required_params = [
        "mp_id_1",
        "mp_id_2",
        "interface_params",
        "functional",
    ]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        print(f"spec is {fw_spec}")
        mpid1 = self.get("mp_id_1")
        mpid2 = self.get("mp_id_2")

        interface_params = self.get("interface_params")
        external_pressure = interface_params.get("external_pressure", 0.0)

        functional = self.get("functional")

        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)

        hl_db = self.get("high_level", True)

        db_high = VaspDB(db_file=db_file, high_level=hl_db)

        slab_surfen_list_1 = fw_spec.get(f"slab_surfen_list_1")
        slab_surfen_list_2 = fw_spec.get(f"slab_surfen_list_2", slab_surfen_list_1)

        # get the surface energies of the slabs
        pairs = list(itertools.product(slab_surfen_list_1, slab_surfen_list_2))
        pairs.sort(key=lambda x: x[0]["surface_energy"] + x[1]["surface_energy"])

        inter_comp_params = get_consolidated_comp_params(
            mpid1=mpid1,
            mpid2=mpid2,
            bulk_coll=f"{functional}.bulk_data",
            db_file=db_file,
            high_level=hl_db,
        )

        for pair in pairs:
            slab1_dict, slab2_dict = pair
            slab1, slab2 = slab1_dict["slab"], slab2_dict["slab"]
            hkl1, hkl2 = slab1_dict["hkl"], slab2_dict["hkl"]
            shift1, shift2 = slab1_dict["shift"], slab2_dict["shift"]

            inter_name = interface_name(mpid1, mpid2, hkl1, hkl2, shift1, shift2)

            inter_data = db_high.find_data(
                collection=functional + ".interface_data",
                fltr={"name": inter_name, "external_pressure": external_pressure},
            )

            if inter_data:
                unrelaxed_structure = inter_data.get("unrelaxed_structure")
                if not unrelaxed_structure:
                    interface_missing = True
                else:
                    interface_missing = False
            else:
                interface_missing = True

            if interface_missing:
                bulk_data_1 = db_high.find_data(
                    collection=f"{functional}.bulk_data", fltr={"mpid": mpid1}
                )
                bulk_data_2 = db_high.find_data(
                    collection=f"{functional}.bulk_data", fltr={"mpid": mpid2}
                )

                bm_1 = bulk_data_1["bulk_moduls"]
                bm_2 = bulk_data_2["bulk_moduls"]

                im = InterfaceMatcher(
                    slab_1=slab1,
                    slab_2=slab2,
                    strain_weight_1=bm_1,
                    strain_weight_2=bm_2,
                    **interface_params,
                )
                top_aligned, bottom_aligned = im.get_centered_slabs()

                if top_aligned and bottom_aligned:
                    interface = im.get_interface()

                    inter_dict = interface.as_dict()
                    bottom_dict = bottom_aligned.as_dict()
                    top_dict = top_aligned.as_dict()

                    db_high.update_data(
                        collection=functional + ".interface_data",
                        fltr={"name": inter_name, "external_pressure": external_pressure},
                        new_values={
                            "$set": {
                                "unrelaxed_structure": inter_dict,
                                "bottom_aligned": bottom_dict,
                                "top_aligned": top_dict,
                                "comp_parameters": inter_comp_params,
                                "interface_params": interface_params,
                            }
                        },
                        upsert=True,
                    )

                    fw_spec["interface_name"] = inter_name
                    return FWAction(update_spec=fw_spec)
                else:
                    # check if we are at the last pair, meaning none of the pairs could be matched
                    if pair == pairs[-1]:
                        return FWAction(defuse_workflow=True)
                    else:
                        continue
            else:
                fw_spec["interface_name"] = inter_name
                return FWAction(update_spec=fw_spec)
