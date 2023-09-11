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
import warnings

from fireworks import (
    explicit_serialize,
    FiretaskBase,
    FWAction,
    Firework,
    Workflow,
)
from pymatgen.core.surface import Slab

from hitmen_utils.db_tools import VaspDB
from surfen.utils.misc_tools import check_input
from surfen.utils.surfen_tools import get_surfen_inputs_from_mpid
from surfen.workflows.subworkflows import surface_energy_swf_from_slab_dict
from triboflow.phys.interface_matcher import InterfaceMatcher
from triboflow.utils.structure_manipulation import flip_slab


@explicit_serialize
class GenerateCandidateSlabs(FiretaskBase):
    """Generate candidate slabs for a given material.

    :param mpid: Material Project's material identifier ID.
    :type mpid: str

    :param material_index: Material index for the given material.
    :type material_index: int

    :param comp_params: Computational parameters for the slab calculations.
    :type comp_params: dict

    :param sg_params: Slab generation parameters.
    :type sg_params: dict

    :param sg_filter: Slab generation filter.
    :type sg_filter: dict

    :param bulk_coll: Name of the bulk collection.
    :type bulk_coll: str

    :param db_file: Path to the database file.
    :type db_file: str

    :param high_level: Name of the high-level database.
    :type high_level: str

    :return: FWAction that updates the spec.
    :rtype: FWAction

    """
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

        nav = VaspDB(db_file=db_file, high_level=high_level)
        bulk_data = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid})
        if bulk_data:
            comp_params_db = bulk_data["comp_parameters"]
            comp_params.update(comp_params_db)

        inputs_list = self.get_inputs_list(
            bulk_coll=bulk_coll,
            comp_params=comp_params,
            db_file=db_file,
            high_level=high_level,
            mpid=mpid,
            sg_filter=sg_filter,
            sg_params=sg_params,
        )
        for input_dict in inputs_list:
            input_dict["mpid"] = mpid
            input_dict["material_index"] = material_index
            input_dict["comp_params"] = comp_params
            input_dict["sg_params"] = sg_params

        fw_spec[f"surfen_inputs_list_{material_index}"] = inputs_list
        return FWAction(update_spec=fw_spec)

    @staticmethod
    def get_inputs_list(
        bulk_coll,
        comp_params,
        db_file,
        high_level,
        mpid,
        sg_filter,
        sg_params,
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
    """Match candidate slabs.

    :param interface_params: Interface parameters.
    :type interface_params: dict

    :return: FWAction that updates the spec.
    :rtype: FWAction

    """
    _fw_name = "Match candidate slabs"
    required_params = []
    optional_params = ["interface_params"]

    def run_task(self, fw_spec):
        interface_params = self.get("interface_params")

        inputs_list_1 = fw_spec["surfen_inputs_list_1"]
        inputs_list_2 = fw_spec.get("surfen_inputs_list_2", inputs_list_1)

        (
            valid_inputs_list_1,
            valid_inputs_list_2,
        ) = self.match_candidate_inputs(
            inputs_list_1, inputs_list_2, interface_params
        )

        if not valid_inputs_list_1 or not valid_inputs_list_2:
            warnings.warn(
                "\n\nNo valid interfaces found for the given slabs. "
                "Defusing the workflow...\n\n"
            )
            return FWAction(defuse_workflow=True)

        fw_spec["valid_surfen_inputs_list_1"] = valid_inputs_list_1
        fw_spec["valid_surfen_inputs_list_2"] = valid_inputs_list_2
        return FWAction(update_spec=fw_spec)

    @staticmethod
    def match_candidate_inputs(
        inputs_list_1, inputs_list_2, interface_params
    ):
        pair_indices = list(
            itertools.product(
                range(len(inputs_list_1)), range(len(inputs_list_2))
            )
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
            interface = im.get_interface()

            if interface is not None:
                valid_pairs.append(pair_index)

        if not valid_pairs:
            return [], []

        # get the unique values of the first index
        unique_indices_1 = list(
            set([pair_index[0] for pair_index in valid_pairs])
        )
        # get the unique values of the second index
        unique_indices_2 = list(
            set([pair_index[1] for pair_index in valid_pairs])
        )
        valid_inputs_list_1 = [
            inputs_list_1[index] for index in unique_indices_1
        ]
        valid_inputs_list_2 = [
            inputs_list_2[index] for index in unique_indices_2
        ]
        return valid_inputs_list_1, valid_inputs_list_2


@explicit_serialize
class RunSurfenSWFOnMatchedSlabs(FiretaskBase):
    """Run the surface energy workflow on the matched slabs.

    :param db_file: Path to the database file.
    :type db_file: str

    :param high_level: Name of the high-level database.
    :type high_level: str

    :param surfen_coll: Name of the surface energy collection.
    :type surfen_coll: str

    :param add_full_relax: Add a full relaxation calculation to the workflow.
    :type add_full_relax: bool

    :return: FWAction that updates the spec.
    """
    _fw_name = (
        "Starts a sub-workflow that calculates surface energies as detour"
    )
    required_params = []
    optional_params = [
        "db_file",
        "high_level",
        "surfen_coll",
        "add_full_relax",
    ]

    def run_task(self, fw_spec):
        db_file = self.get("db_file")
        high_level = self.get("high_level")
        surfen_coll = self.get("surfen_coll")
        add_full_relax = self.get("add_full_relax")

        inputs_list_1 = fw_spec[f"valid_surfen_inputs_list_1"]
        inputs_list_2 = fw_spec[f"valid_surfen_inputs_list_2"]

        inputs_list = inputs_list_1 + inputs_list_2
        unique_inputs_list = list(
            {inp["uid"]: inp for inp in inputs_list}.values()
        )

        fw_spec["surfen_inputs_combined_list"] = unique_inputs_list

        wfs = []
        for slab_dict in unique_inputs_list:
            mpid = slab_dict["mpid"]
            comp_params = slab_dict["comp_params"]
            sg_params = slab_dict["sg_params"]
            wf = surface_energy_swf_from_slab_dict(
                mpid=mpid,
                slab_dict=slab_dict,
                surfen_coll=surfen_coll,
                db_file=db_file,
                high_level=high_level,
                sg_params=sg_params,
                comp_params=comp_params,
                add_full_relax=add_full_relax,
            )
            if wf:
                wfs.append(wf)
        return FWAction(detours=wfs, update_spec=fw_spec)


@explicit_serialize
class GetSlabSurfenListFromUids(FiretaskBase):
    """Query the surface energies with the list of uids provided.

    :param db_file: Path to the database file.
    :type db_file: str

    :param high_level: Name of the high-level database.
    :type high_level: str

    :param surfen_coll: Name of the surface energy collection.
    :type surfen_coll: str

    :return: FWAction that updates the spec.
    """
    _fw_name = "Query the surface energies with the list of uids provided"
    required_params = []
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

        surfen_inputs_combined_list = fw_spec["surfen_inputs_combined_list"]
        surfen_inp_dict = {
            inp["uid"]: inp for inp in surfen_inputs_combined_list
        }
        uids = [inp["uid"] for inp in surfen_inputs_combined_list]
        nav = VaspDB(db_file=inp["db_file"], high_level=inp["high_level"])
        results = list(
            nav.find_many_data(
                inp["surfen_coll"], {"uid": {"$in": uids}}, {"calcs": 0}
            )
        )

        slab_surfen_list = []
        for result in results:
            uid = result["uid"]
            slab_dict = surfen_inp_dict[uid]
            mpid = slab_dict["mpid"]
            material_index = slab_dict["material_index"]
            slab = Slab.from_dict(result["structure"])
            shift = round(slab.shift, 2)
            surface_energies = result["surface_energy"]
            slab_params = result["slab_params"]
            hkl = result["hkl"]
            terminations = result["terminations"]
            surfen_min = min(surface_energies.values())

            if not slab_params.get("sym", False):
                surfen_min_key = min(
                    surface_energies, key=surface_energies.get
                )
                if surfen_min_key == "bottom":
                    slab = flip_slab(slab)
                    termination_min = terminations["bottom"]
                    # flip the terminations
                    terminations["bottom"], terminations["top"] = (
                        terminations["top"],
                        terminations["bottom"],
                    )
                    # flip the surface energies
                    surface_energies["bottom"], surface_energies["top"] = (
                        surface_energies["top"],
                        surface_energies["bottom"],
                    )
                else:
                    termination_min = terminations["top"]
            else:
                termination_min = terminations["top"]

            slab_surfen_list.append(
                {
                    "mpid": mpid,
                    "material_index": material_index,
                    "slab": slab,
                    "hkl": hkl,
                    "shift": shift,
                    "termination": termination_min,
                    "surface_energy": surfen_min,
                    "slab_params": slab_params,
                    "surface_energies": surface_energies,
                    "terminations": terminations,
                }
            )

        slab_surfen_list.sort(key=lambda x: x["surface_energy"])

        slab_surfen_list_1 = [
            entry for entry in slab_surfen_list if entry["material_index"] == 1
        ]

        slab_surfen_list_2 = [
            entry for entry in slab_surfen_list if entry["material_index"] == 2
        ]

        # if either of them is an empty list, set it to the other one
        # TODO: FIND A BETTER WAY TO DO THIS, PERHAPS IN MakeHeterostructure
        if not slab_surfen_list_1:
            slab_surfen_list_1 = slab_surfen_list_2
        if not slab_surfen_list_2:
            slab_surfen_list_2 = slab_surfen_list_1
        fw_spec[f"slab_surfen_list_1"] = slab_surfen_list_1
        fw_spec[f"slab_surfen_list_2"] = slab_surfen_list_2

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class GetCandidatesForHeteroStructure(FiretaskBase):
    _fw_name = "Match Slabs and Calculate Surface Energies"
    required_params = [
        "mpid_1",
        "mpid_2",
        "comp_params_1",
        "comp_params_2",
        "interface_params",
        "sg_params_1",
        "sg_params_2",
        "sg_filter_1",
        "sg_filter_2",
    ]
    optional_params = [
        "db_file",
        "high_level",
        "surfen_coll",
        "bulk_coll",
        "add_full_relax",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}

        mpid_1 = inp["mpid_1"]
        mpid_2 = inp["mpid_2"]
        comp_params_1 = inp["comp_params_1"]
        comp_params_2 = inp["comp_params_2"]
        interface_params = inp["interface_params"]
        sg_params_1 = inp["sg_params_1"]
        sg_params_2 = inp["sg_params_2"]
        sg_filter_1 = inp["sg_filter_1"]
        sg_filter_2 = inp["sg_filter_2"]
        db_file = inp["db_file"]
        high_level = inp["high_level"]
        surfen_coll = inp["surfen_coll"]
        bulk_coll = inp["bulk_coll"]
        add_full_relax = inp["add_full_relax"]

        fw_list = []

        fw1 = Firework(
            [
                GenerateCandidateSlabs(
                    mpid=mpid_1,
                    material_index=1,
                    comp_params=comp_params_1,
                    sg_params=sg_params_1,
                    sg_filter=sg_filter_1,
                    bulk_coll=bulk_coll,
                    db_file=db_file,
                    high_level=high_level,
                )
            ],
            name=f"Generate Candidate Slabs for {mpid_1}",
        )

        fw_list.append(fw1)

        if (
            (mpid_2 != mpid_1)
            or (sg_params_1 != sg_params_2)
            or (sg_filter_1 != sg_filter_2)
        ):
            fw2 = Firework(
                [
                    GenerateCandidateSlabs(
                        mpid=mpid_2,
                        material_index=2,
                        comp_params=comp_params_2,
                        sg_params=sg_params_2,
                        sg_filter=sg_filter_2,
                        bulk_coll=bulk_coll,
                        db_file=db_file,
                        high_level=high_level,
                    )
                ],
                name=f"Generate Candidate Slabs for {mpid_2}",
            )
            fw_list.append(fw2)
        else:
            fw2 = None

        parents = [fw1, fw2] if fw2 else [fw1]

        fw3 = Firework(
            [MatchCandidateSlabs(interface_params=interface_params)],
            name="Match Candidate Slabs",
            parents=parents,
        )

        fw_list.append(fw3)

        fw4 = Firework(
            [
                RunSurfenSWFOnMatchedSlabs(
                    db_file=db_file,
                    high_level=high_level,
                    surfen_coll=surfen_coll,
                    add_full_relax=add_full_relax,
                )
            ],
            name="Start Surface Energy from Slab Dicts",
            parents=[fw3],
        )

        fw_list.append(fw4)

        fw5 = Firework(
            [
                GetSlabSurfenListFromUids(
                    db_file=db_file,
                    high_level=high_level,
                    surfen_coll=surfen_coll,
                )
            ],
            name="Get Slab Surfen List from Uids",
            parents=[fw4],
        )

        fw_list.append(fw5)

        wf = Workflow(
            fw_list, name=f"Match Slabs and Calculate Surface Energies"
        )
        return FWAction(detours=wf, update_spec=fw_spec)


# @explicit_serialize
# class RunSurfenSwfGetEnergies(FiretaskBase):
#     _fw_name = "Run the surface energy workflow and get the energies"
#     required_params = ["mpid", "material_index"]
#     optional_params = [
#         "comp_params",
#         "sg_params",
#         "sg_filter",
#         "db_file",
#         "high_level",
#         "custom_id",
#         "surfen_coll",
#         "bulk_coll",
#         "add_full_relax",
#     ]
#
#     def run_task(self, fw_spec):
#         req_inp = {k: self.get(k) for k in self.required_params}
#         opt_inp = {k: self.get(k) for k in self.optional_params}
#         opt_inp = check_input(opt_inp, self.optional_params)
#         inp = {**req_inp, **opt_inp}
#         mpid = inp["mpid"]
#         comp_params = inp["comp_params"]
#         sg_params = inp["sg_params"]
#         db_file = inp["db_file"]
#         high_level = inp["high_level"]
#         # create a name_suffix variable and use miller or max_index from sg_params, whichever is available
#         name_suffix = ""
#         if sg_params.get("miller", False):
#             name_suffix = (
#                 f"Miller {sg_params['miller']} - Material Index {inp['material_index']}"
#             )
#         elif sg_params.get("max_index", False):
#             name_suffix = (
#                 f"MMI {sg_params['max_index']} - Material Index {inp['material_index']}"
#             )
#         # check if comp_params_loc is given
#         nav = VaspDB(db_file=db_file, high_level=high_level)
#         bulk_coll = inp["bulk_coll"]
#
#         bulk_data = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid})
#         if bulk_data:
#             comp_params_db = bulk_data["comp_parameters"]
#             comp_params.update(comp_params_db)
#
#         fw1 = Firework(
#             [
#                 StartSurfaceEnergy(
#                     mpid=mpid,
#                     comp_params=comp_params,
#                     sg_params=inp["sg_params"],
#                     sg_filter=inp["sg_filter"],
#                     db_file=db_file,
#                     high_level=high_level,
#                     custom_id=inp["custom_id"],
#                     surfen_coll=inp["surfen_coll"],
#                     bulk_coll=inp["bulk_coll"],
#                     add_full_relax=True,
#                 )
#             ],
#             name=f"Run Surface Energy Workflow for {mpid} {name_suffix}",
#         )
#
#         fw2 = Firework(
#             [
#                 GetSlabSurfenListFromUids(
#                     material_index=inp["material_index"],
#                     db_file=db_file,
#                     high_level=high_level,
#                     surfen_coll=inp["surfen_coll"],
#                 )
#             ],
#             name=f"Get slabs and corresponding surface energies for {mpid} {name_suffix}",
#             parents=fw1,
#         )
#
#         wf = Workflow(
#             [fw1, fw2],
#             name=f"Run Surface Energy Workflow and Extract Surface Energies for {mpid} {name_suffix}",
#         )
#         return FWAction(detours=wf, update_spec=fw_spec)
