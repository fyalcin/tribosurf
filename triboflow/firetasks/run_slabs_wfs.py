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

import itertools

from fireworks import (
    explicit_serialize,
    FiretaskBase,
    FWAction,
    Workflow,
    Firework,
)
from pymatgen.core.surface import Slab

from hitmen_utils.db_tools import VaspDB
from surfen.firetasks.start_swfs import StartSurfaceEnergy
from surfen.utils.misc_tools import check_input
from surfen.utils.surfen_tools import get_surfen_inputs_from_mpid
from surfen.workflows.subworkflows import surface_energy_swf_from_slab_dict
from triboflow.phys.interface_matcher import InterfaceMatcher
from triboflow.utils.structure_manipulation import flip_slab


@explicit_serialize
class GenerateCandidateSlabs(FiretaskBase):
    _fw_name = "Generate candidate slabs"
    required_params = ["mpid", "material_index"]
    optional_params = [
        "comp_params",
        "sg_params",
        "sg_filter",
        "bulk_coll",
        "db_file",
        "high_level",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        mpid = inp["mpid"]
        sg_params = inp["sg_params"]
        sg_filter = inp["sg_filter"]
        bulk_coll = inp["bulk_coll"]
        db_file = inp["db_file"]
        high_level = inp["high_level"]
        material_index = inp["material_index"]
        comp_params = inp["comp_params"]

        inputs_list = self.get_inputs_list(
            bulk_coll=bulk_coll,
            comp_params=comp_params,
            db_file=db_file,
            high_level=high_level,
            mpid=mpid,
            sg_filter=sg_filter,
            sg_params=sg_params,
        )

        fw_spec[f"surfen_inputs_list_{material_index}"] = inputs_list
        return FWAction(update_spec=fw_spec)

    @staticmethod
    def get_inputs_list(
        bulk_coll, comp_params, db_file, high_level, mpid, sg_filter, sg_params
    ):
        inputs_list = get_surfen_inputs_from_mpid(
            mpid=mpid,
            bulk_coll=bulk_coll,
            sg_params=sg_params,
            sg_filter=sg_filter,
            comp_params=comp_params,
            db_file=db_file,
            high_level=high_level,
            custom_id=None,
        )
        return inputs_list


@explicit_serialize
class MatchCandidateSlabs(FiretaskBase):
    _fw_name = "Match candidate slabs"
    required_params = []
    optional_params = ["interface_params"]

    def run_task(self, fw_spec):
        interface_params = self.get("interface_params")

        inputs_list_1 = fw_spec["surfen_inputs_list_1"]
        inputs_list_2 = fw_spec.get("surfen_inputs_list_2", inputs_list_1)

        unique_inputs_list_1, unique_inputs_list_2 = self.match_candidate_inputs(
            inputs_list_1, inputs_list_2, interface_params
        )

        fw_spec["surfen_inputs_list_1"] = unique_inputs_list_1
        fw_spec["surfen_inputs_list_2"] = unique_inputs_list_2

        return FWAction(update_spec=fw_spec)

    @staticmethod
    def match_candidate_inputs(inputs_list_1, inputs_list_2, interface_params):
        pair_indices = list(
            itertools.product(range(len(inputs_list_1)), range(len(inputs_list_2)))
        )
        valid_pairs = []
        for pair_index in pair_indices:
            slab_dict_1 = inputs_list_1[pair_index[0]]
            slab_dict_2 = inputs_list_2[pair_index[1]]

            slab_1 = slab_dict_1["struct"]
            slab_2 = slab_dict_2["struct"]

            im = InterfaceMatcher(
                slab_1=slab_1,
                slab_2=slab_2,
                strain_weight_1=1,
                strain_weight_2=1,
                **interface_params,
            )
            top_aligned, bottom_aligned = im.get_centered_slabs()

            if top_aligned and bottom_aligned:
                valid_pairs.append(pair_index)
        # get the unique values of the first index
        unique_indices_1 = list(set([pair_index[0] for pair_index in valid_pairs]))
        # get the unique values of the second index
        unique_indices_2 = list(set([pair_index[1] for pair_index in valid_pairs]))
        valid_inputs_list_1 = [inputs_list_1[index] for index in unique_indices_1]
        valid_inputs_list_2 = [inputs_list_2[index] for index in unique_indices_2]
        return valid_inputs_list_1, valid_inputs_list_2


@explicit_serialize
class GetSlabSurfenListFromUids(FiretaskBase):
    _fw_name = "Query the surface energies with the list of uids provided"
    required_params = ["material_index"]
    optional_params = [
        "db_file",
        "high_level",
        "surfen_coll",
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
    required_params = ["mpid", "material_index"]
    optional_params = [
        "comp_params",
        "sg_params",
        "sg_filter",
        "db_file",
        "high_level",
        "custom_id",
        "surfen_coll",
        "bulk_coll",
        "add_full_relax",
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
        nav = VaspDB(db_file=db_file, high_level=high_level)
        bulk_coll = inp["bulk_coll"]

        bulk_data = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid})
        if bulk_data:
            comp_params_db = bulk_data["comp_parameters"]
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
                    material_index=inp["material_index"],
                    db_file=db_file,
                    high_level=high_level,
                    surfen_coll=inp["surfen_coll"],
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


class MatchSlabsCalculateSurfen(FiretaskBase):
    _fw_name = "Run the surface energy workflow on matched slabs"
    required_params = ["mpid_1", "mpid_2"]
    optional_params = [
        "comp_params",
        "interface_params",
        "sg_params_1",
        "sg_params_2",
        "sg_filter_1",
        "sg_filter_2",
        "db_file",
        "high_level",
        "custom_id",
        "surfen_coll",
        "bulk_coll",
        "add_full_relax",
    ]

    def run_task(self, fw_spec):
        mpid_1 = self.get("mpid_1")
        mpid_2 = self.get("mpid_2")
        comp_params = self.get("comp_params")
        interface_params = self.get("interface_params")
        sg_params_1 = self.get("sg_params_1")
        sg_params_2 = self.get("sg_params_2")
        sg_filter_1 = self.get("sg_filter_1")
        sg_filter_2 = self.get("sg_filter_2")
        db_file = self.get("db_file")
        high_level = self.get("high_level")
        custom_id = self.get("custom_id")
        surfen_coll = self.get("surfen_coll")
        bulk_coll = self.get("bulk_coll")
        add_full_relax = self.get("add_full_relax")

        wf_list = []

        fw1 = Firework(GenerateCandidateSlabs(
            mpid=mpid_1,
            material_index=1,
            comp_params=comp_params,
            sg_params=sg_params_1,
            sg_filter=sg_filter_1,
            bulk_coll=bulk_coll,
            db_file=db_file,
            high_level=high_level,
        ),
            name=f"Generate candidate slabs for {mpid_1}",
        )
        wf_list.append(fw1)

        if (
                (mpid_2 != mpid_1)
                or (sg_params_1 != sg_params_2)
                or (sg_filter_1 != sg_filter_2)
        ):
            fw2 = Firework(GenerateCandidateSlabs(
                mpid=mpid_2,
                material_index=2,
                comp_params=comp_params,
                sg_params=sg_params_2,
                sg_filter=sg_filter_2,
                bulk_coll=bulk_coll,
                db_file=db_file,
                high_level=high_level,
            ),
                name=f"Generate candidate slabs for {mpid_2}",
            )
            wf_list.append(fw2)
        else:
            fw2 = None

        fw3 = Firework(MatchCandidateSlabs(
            interface_params=interface_params),
            name=f"Match candidate slabs for {mpid_1} and {mpid_2}",
            parents=[fw1, fw2] if fw2 else [fw1],
        )

        wf_list.append(fw3)


@explicit_serialize
class StartSurfaceEnergyFromSlabDicts(FiretaskBase):
    _fw_name = "Starts a sub-workflow that calculates surface energies as detour"
    required_params = ["mpid", "material_index"]
    optional_params = [
        "comp_params",
        "sg_params",
        "db_file",
        "high_level",
        "surfen_coll",
        "add_full_relax",
    ]

    def run_task(self, fw_spec):
        mpid = self.get("mpid")
        material_index = self.get("material_index")
        comp_params = self.get("comp_params")
        sg_params = self.get("sg_params")
        db_file = self.get("db_file")
        high_level = self.get("high_level")
        surfen_coll = self.get("surfen_coll")
        add_full_relax = self.get("add_full_relax")

        inputs_list = fw_spec[f"surfen_inputs_list_{material_index}"]
        wfs = []
        for slab_dict in inputs_list:
            wf = surface_energy_swf_from_slab_dict(
                mpid=mpid,
                slab_dict=slab_dict,
                surfen_coll=surfen_coll,
                db_file=db_file,
                high_level=high_level,
                sg_params=sg_params,
                comp_params=comp_params,
                add_full_relax=add_full_relax)
        return FWAction(detours=wfs, update_spec=fw_spec)

