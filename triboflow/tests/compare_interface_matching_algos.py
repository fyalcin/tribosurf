#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:42:27 2022

@author: wolloch
"""
from itertools import combinations

from pymatgen.core.surface import Slab

from triboflow.phys.interface_matcher import InterfaceMatcher
from triboflow.utils.database import Navigator
from triboflow.utils.structure_manipulation import interface_name

interface_parameters = {
    "max_area": 100,
    "interface_distance": "auto",
    "max_angle_diff": 1,
    "max_mismatch": 0.05,
    "r1r2_tol": 0.5,
    "best_match": "area",
    "vacuum": 15,
}

nav = Navigator(high_level="production")
slab_list = []
for slab_dict in nav.find_many_data("PBE.slab_data", {}):
    try:
        slab_list.append(
            {
                "mpid": slab_dict["mpid"],
                "miller": slab_dict["miller"],
                "slab": Slab.from_dict(slab_dict["relaxed_slab"]),
            }
        )
    except:
        pass

results = {
    "statistics": {
        "nr_of_interfaces": 0,
        "failed_matches": 0,
        "both_methods_matched": 0,
        "both_methods_matched_differently": 0,
        "both_methods_matched_same": 0,
        "pmg_only_matched": 0,
        "mpi_only_matched": 0,
    },
    "parameters": interface_parameters,
    "matched_differently": [],
    "pmg_only": [],
    "mpi_only": [],
}

slab_comb_list = list(combinations(slab_list, 2))

for i, slab_combinations in enumerate(slab_comb_list):
    progress = int((i + 1) / len(slab_comb_list) * 100)
    results["statistics"]["nr_of_interfaces"] += 1

    name = interface_name(
        slab_combinations[0]["mpid"],
        slab_combinations[0]["miller"],
        slab_combinations[1]["mpid"],
        slab_combinations[1]["miller"],
    )
    print("-".center(100, "-"))
    print(f"{i+1} of {len(slab_comb_list)}".center(progress, "|"))
    print("-".center(100, "-"))
    IM_pmg = InterfaceMatcher(
        slab_combinations[0]["slab"],
        slab_combinations[1]["slab"],
        implementation="pymatgen",
        **interface_parameters,
    )
    IM_mpi = InterfaceMatcher(
        slab_combinations[0]["slab"],
        slab_combinations[1]["slab"],
        implementation="MPInterfaces",
        **interface_parameters,
    )
    interface_pmg = IM_pmg.get_interface()
    interface_mpi = IM_mpi.get_interface()
    if interface_pmg and interface_mpi:
        results["statistics"]["both_methods_matched"] += 1
        if interface_pmg == interface_mpi:
            results["statistics"]["both_methods_matched_same"] += 1
        else:
            results["statistics"]["both_methods_matched_differently"] += 1
            results["matched_differently"].append(name)
            interface_pmg.to("poscar", "./test_interfaces/" + name + "_pmg.vasp")
            interface_mpi.to("poscar", "./test_interfaces/" + name + "_mpi.vasp")
        results[name] = {
            "pmg": {
                "interface": interface_pmg,
                "area": interface_pmg.interface_properties["area"],
                "strain": interface_pmg.interface_properties["strain"],
                "lattice_abc": interface_pmg.lattice.abc,
            },
            "mpi": {
                "interface": interface_mpi,
                "area": interface_mpi.interface_properties["area"],
                "strain": interface_mpi.interface_properties["strain"],
                "lattice_abc": interface_mpi.lattice.abc,
            },
        }

    elif interface_pmg:
        results["statistics"]["pmg_only_matched"] += 1
        results["pmg_only"].append(name)
        results[name] = {
            "pmg": {
                "interface": interface_pmg,
                "area": interface_pmg.interface_properties["area"],
                "strain": interface_pmg.interface_properties["strain"],
                "lattice_abc": interface_pmg.lattice.abc,
            },
            "mpi": None,
        }
        interface_pmg.to("poscar", "./test_interfaces/" + name + "_pmg.vasp")
    elif interface_mpi:
        results["statistics"]["mpi_only_matched"] += 1
        results["mpi_only"].append(name)
        results[name] = {
            "mpi": {
                "interface": interface_mpi,
                "area": interface_mpi.interface_properties["area"],
                "strain": interface_mpi.interface_properties["strain"],
                "lattice_abc": interface_mpi.lattice.abc,
            },
            "pmg": None,
        }
        interface_mpi.to("poscar", "./test_interfaces/" + name + "_mpi.vasp")
    else:
        results["statistics"]["failed_matches"] += 1
        results[name] = {"pmg": None, "mpi": None}

stats = results["statistics"]
results["statistics"]["percentage_matched_both"] = (
    stats["both_methods_matched"] / stats["nr_of_interfaces"] * 100
)
results["statistics"]["percentage_matched_pmg"] = (
    (stats["both_methods_matched"] + stats["pmg_only_matched"])
    / stats["nr_of_interfaces"]
    * 100
)
results["statistics"]["percentage_matched_mpi"] = (
    (stats["both_methods_matched"] + stats["mpi_only_matched"])
    / stats["nr_of_interfaces"]
    * 100
)
