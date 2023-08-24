# -*- coding: utf-8 -*-
import pprint
import argparse
from pymatgen.io.vasp.inputs import Poscar
from triboflow.phys.interface_matcher import InterfaceMatcher
from triboflow.utils.structure_manipulation import slab_from_file

"""
Match interfaces using Slabs from POSCAR files.
"""


def match_the_interface(slab_1, slab_2, inter_params={}):
    max_area = inter_params.get("max_area", 500)
    max_length_tol = inter_params.get("max_length_tol", 0.01)
    max_angle_tol = inter_params.get("max_angle_tol", 0.01)
    max_area_ratio_tol = inter_params.get("max_area_ratio_tol", 0.1)
    separation = inter_params.get("separation", "auto")

    IM = InterfaceMatcher(
        slab_1,
        slab_2,
        max_area=max_area,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
        max_area_ratio_tol=max_area_ratio_tol,
        interface_distance=separation,
    )

    top_aligned, bottom_aligned = IM.get_centered_slabs()

    if top_aligned and bottom_aligned:
        interface = IM.get_interface()
        return_dict = {
            "interface": interface,
            "matched_bottom_slab": bottom_aligned,
            "matched_top_slab": top_aligned,
        }
        return return_dict
    else:
        print(
            "No match of the given slabs could be found with the "
            "given parameters:"
        )
        pprint.pprint(inter_params)
        print("")
        return None


def GetUserInput():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-bs",
        "--bottomslab",
        dest="slab_name_1",
        type=str,
        help="Filename of the bottom slab POSCAR file",
    )
    parser.add_argument(
        "-bm",
        "--bottom_miller",
        dest="miller_1",
        nargs=3,
        type=int,
        help="3 Miller indices of bottom slab surface",
    )
    parser.add_argument(
        "-ts",
        "--topslab",
        dest="slab_name_2",
        type=str,
        help="Filename of the top slab POSCAR file",
    )
    parser.add_argument(
        "-tm",
        "--top_miller",
        dest="miller_2",
        nargs=3,
        type=int,
        help="3 Miller indices of top slab surface",
    )
    parser.add_argument(
        "-marea",
        "--max_area",
        dest="max_area",
        type=float,
        default=200.0,
        help="Maximally allowed cell cross section area",
    )
    parser.add_argument(
        "-mat",
        "--max_angle_tol",
        dest="max_angle_tol",
        type=float,
        default=0.01,
        help="Maximal relative missmatch for angles",
    )
    parser.add_argument(
        "-mlt",
        "--max_length_tol",
        dest="max_length_tol",
        type=float,
        default=0.01,
        help="Maximal relative missmatch in vector length",
    )
    parser.add_argument(
        "-rt",
        "--max_area_ratio_tol",
        dest="max_area_ratio_tol",
        type=float,
        default=0.1,
        help="Maximal relative difference in area ratio.",
    )
    parser.add_argument(
        "-s",
        "--separation",
        dest="separation",
        type=float,
        default=2.5,
        help="Distance between slabs in Angstrom.",
    )
    args = parser.parse_args()
    return args


def WriteOutput(interface_dict):
    bottom_slab = interface_dict["matched_bottom_slab"]
    top_slab = interface_dict["matched_top_slab"]
    interface = interface_dict["interface"]
    Poscar(bottom_slab).write_file("Matched_bottom_slab.vasp")
    Poscar(top_slab).write_file("Matched_top_slab.vasp")
    Poscar(interface).write_file("Matched_interface.vasp")
    return


if __name__ == "__main__":
    args = GetUserInput()
    slab_1 = slab_from_file(args.miller_1, args.slab_name_1)
    slab_2 = slab_from_file(args.miller_2, args.slab_name_2)
    inter_params = {
        "max_area": args.max_area,
        "max_angle_diff": args.max_angle_diff,
        "max_missmatch": args.max_missmatch,
        "r1r2_tol": args.r1r2_tol,
        "separation": args.separation,
    }
    Matching_info = match_the_interface(slab_1, slab_2, inter_params)
    if Matching_info:
        WriteOutput(Matching_info)
