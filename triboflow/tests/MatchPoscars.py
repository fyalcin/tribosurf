# -*- coding: utf-8 -*-
import pprint
import argparse
from pymatgen.io.vasp.inputs import Poscar
from mpinterfaces.utils import slab_from_file
from mpinterfaces.transformations import get_aligned_lattices, \
    get_interface
"""
Match interfaces using Slabs from POSCAR files.
"""


def MatchTheInterface(slab_1, slab_2, inter_params={}):
    
    max_area = inter_params.get('max_area', 500)
    max_missmatch = inter_params.get('max_missmatch', 0.01)
    max_angle_diff = inter_params.get('max_angle_diff', 2.0)
    r1r2_tol = inter_params.get('r1r2_tol', 0.05)
    separation = inter_params.get('separation', 2.5)
    
    bottom_aligned, top_aligned = get_aligned_lattices(
                                        slab_1,
                                        slab_2,
                                        max_area = max_area,
                                        max_mismatch = max_missmatch,
                                        max_angle_diff = max_angle_diff,
                                        r1r2_tol = r1r2_tol)
    if bottom_aligned:
        hetero_interfaces = get_interface(bottom_aligned,
                                      top_aligned,
                                      nlayers_2d = 1,
                                      nlayers_substrate = 1,
                                      separation = separation)
        return_dict = {"interface": hetero_interfaces,
                       "matched_bottom_slab": bottom_aligned,
                       "matched_top_slab": top_aligned}
        return return_dict
    else:
        print('No match of the given slabs could be found with the '
              'given parameters:')
        pprint.pprint(inter_params)
        print('')
        return None
    
def GetUserInput():
    parser = argparse.ArgumentParser()
    parser.add_argument("-bs", "--bottomslab", dest = "slab_name_1",
                        type = str,
                        help="Filename of the bottom slab POSCAR file")
    parser.add_argument("-bm", "--bottom_miller", dest = "miller_1",
                        nargs = 3, type = int,
                        help="3 Miller indices of bottom slab surface")
    parser.add_argument("-ts", "--topslab", dest = "slab_name_2",
                        type = str,
                        help="Filename of the top slab POSCAR file")
    parser.add_argument("-tm", "--top_miller", dest = "miller_2",
                        nargs = 3, type = int,
                        help="3 Miller indices of top slab surface")
    parser.add_argument("-marea", "--max_area", dest = "max_area",
                        type = float, default = 200.0,
                        help = "Maximally allowed cell cross section area")
    parser.add_argument("-mangle", "--max_angle_diff", dest = "max_angle_diff",
                        type = float, default = 1.0,
                        help = "Maximally allowed missmatch angle difference")
    parser.add_argument("-mm", "--max_missmatch", dest = "max_missmatch",
                        type = float, default = 0.05,
                        help = "Maximally allowed missmatch")
    parser.add_argument("-rt", "--r1r2_tol", dest = "r1r2_tol",
                        type = float, default = 0.2,
                        help = "Tolerance parameter for matching lattices.")
    parser.add_argument("-s", "--separation", dest = "separation",
                        type = float, default = 2.5,
                        help = "Distance between slabs in Angstrom.")
    args = parser.parse_args()
    return args

def WriteOutput(interface_dict):
    bottom_slab = interface_dict['matched_bottom_slab']
    top_slab = interface_dict['matched_top_slab']
    interface = interface_dict['interface']
    Poscar(bottom_slab).write_file('Matched_bottom_slab.vasp')
    Poscar(top_slab).write_file('Matched_top_slab.vasp')
    Poscar(interface).write_file('Matched_interface.vasp')
    return

if __name__ == "__main__":
    args = GetUserInput()
    slab_1 = slab_from_file(args.miller_1, args.slab_name_1)
    slab_2 = slab_from_file(args.miller_2, args.slab_name_2)
    inter_params = {'max_area': args.max_area,
                    'max_angle_diff': args.max_angle_diff,
                    'max_missmatch': args.max_missmatch,
                    'r1r2_tol': args.r1r2_tol,
                    'separation': args.separation}
    Matching_info = MatchTheInterface(slab_1, slab_2, inter_params)
    if Matching_info:
        WriteOutput(Matching_info)
    