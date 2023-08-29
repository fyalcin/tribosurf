#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# TODO - change copyright
"""

Collection of Firetasks to start workflows concerning surfaces.

The module contains the following Firetasks:

    - GetSlabSurfenListFromUids: query the surface energies with the list of uids provided
    - RunSurfenSwfGetEnergies: run the surface energy workflow and get the surface energies

"""

__author__ = "Firat Yalcin"
__copyright__ = ""
__contact__ = "firat.yalcin@univie.ac.at"

from fireworks import (
    explicit_serialize,
    FiretaskBase,
    FWAction,
    Workflow,
    Firework,
)
from pymatgen.core.surface import Slab

from hitmen_utils.db_tools import VaspDB
from surfen.utils.misc_tools import check_input
from surfen.firetasks.start_swfs import StartSurfaceEnergy
from triboflow.utils.structure_manipulation import flip_slab


@explicit_serialize
class GetSlabSurfenListFromUids(FiretaskBase):
    _fw_name = "Query the surface energies with the list of uids provided"
    required_params = []
    optional_params = [
        "db_file",
        "high_level",
        "fake_calc",
        "surfen_coll",
        "material_index",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        print(f"inp: {inp}")
        prev_swf_info = fw_spec["StartSurfaceEnergy_info"]
        mpid, uids = prev_swf_info["mpid"], prev_swf_info["uids"]
        nav = VaspDB(db_file=inp["db_file"], high_level=inp["high_level"])
        results = list(
            nav.find_many_data(inp["surfen_coll"], {"uid": {"$in": uids}}, {"calcs": 0})
        )
        # for result in results:
        #     slab = Slab.from_dict(result["structure"])
        #     hkl = result["hkl"]
        #     uid = result["uid"]
        #     nav.update_data(
        #         collection="test",
        #         fltr={"mpid": mpid},
        #         new_values={
        #             "$set": {
        #                 uid + ".hkl": hkl,
        #                 uid + ".shift": slab.shift,
        #                 uid + ".surface_energy": result["surface_energy"],
        #                 uid + ".structure": result["structure"],
        #                 uid + ".terminations": result["terminations"],
        #             }
        #         },
        #         upsert=True,
        #     )

        slab_surfen_list = []
        for slab_dict in results:
            slab = Slab.from_dict(slab_dict["structure"])
            shift = round(slab.shift, 2)
            surfen_dict = slab_dict["surface_energy"]
            slab_params = slab_dict["slab_params"]
            hkl = slab_dict["hkl"]
            terminations = slab_dict["terminations"]
            surfen_min = min(surfen_dict.values())

            if not slab_params.get("sym", False):
                surfen_min_key = min(surfen_dict, key=surfen_dict.get)
                if surfen_min_key == "bottom":
                    slab = flip_slab(slab)
                    termination_min = terminations["bottom"]
                else:
                    termination_min = terminations["top"]
            else:
                termination_min = terminations["top"]

            slab_surfen_list.append(
                {
                    "slab": slab,
                    "hkl": hkl,
                    "shift": shift,
                    "termination": termination_min,
                    "surface_energy": surfen_min,
                    "slab_params": slab_params,
                }
            )

        slab_surfen_list.sort(key=lambda x: x["surface_energy"])
        fw_spec[f"slab_surfen_list_{inp['material_index']}"] = slab_surfen_list

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class RunSurfenSwfGetEnergies(FiretaskBase):
    _fw_name = "Run the surface energy workflow and get the energies"
    required_params = ["mpid"]
    optional_params = [
        "comp_params",
        "comp_params_from_db",
        "sg_params",
        "sg_filter",
        "db_file",
        "high_level",
        "custom_id",
        "surfen_coll",
        "bulk_coll",
        "add_full_relax",
        "material_index",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        mpid = inp["mpid"]
        comp_params = inp["comp_params"]
        sg_params = inp["sg_params"]
        db_file = inp["db_file"]
        high_level = inp["high_level"]
        # create a name_suffix variable and use miller or max_index from sg_params, whichever is available
        name_suffix = ""
        if sg_params.get("miller", False):
            name_suffix = (
                f"Miller {sg_params['miller']} - Material Index {inp['material_index']}"
            )
        elif sg_params.get("max_index", False):
            name_suffix = (
                f"MMI {sg_params['max_index']} - Material Index {inp['material_index']}"
            )
        # check if comp_params_loc is given
        if inp["comp_params_from_db"]:
            nav = VaspDB(db_file=db_file, high_level=high_level)
            bulk_coll = inp["bulk_coll"]

            bulk_data = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid})
            if bulk_data:
                comp_params_db = bulk_data["comp_parameters"]
            else:
                comp_params_db = None

            if comp_params_db:
                comp_params.update(comp_params_db)

        fw1 = Firework(
            [
                StartSurfaceEnergy(
                    mpid=mpid,
                    comp_params=comp_params,
                    sg_params=inp["sg_params"],
                    sg_filter=inp["sg_filter"],
                    db_file=db_file,
                    high_level=high_level,
                    custom_id=inp["custom_id"],
                    surfen_coll=inp["surfen_coll"],
                    bulk_coll=inp["bulk_coll"],
                    add_full_relax=True,
                )
            ],
            name=f"Run Surface Energy Workflow for {mpid} {name_suffix}",
        )

        fw2 = Firework(
            [
                GetSlabSurfenListFromUids(
                    db_file=db_file,
                    high_level=high_level,
                    surfen_coll=inp["surfen_coll"],
                    material_index=inp["material_index"],
                )
            ],
            name=f"Get slabs and corresponding surface energies for {mpid} {name_suffix}",
            parents=fw1,
        )

        wf = Workflow(
            [fw1, fw2],
            name=f"Run Surface Energy Workflow and Extract Surface Energies for {mpid} {name_suffix}",
        )
        return FWAction(detours=wf, update_spec=fw_spec)
