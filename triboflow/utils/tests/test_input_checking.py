import unittest
from copy import deepcopy
from triboflow.utils.check_WF_inputs import check_hetero_wf_inputs

inputs_correct = {
    "material": {
        "formula": "Si",
        "mpid": "mp-149",
    },
    "comp_params": {
        "functional": "PBE",
        "volume_tolerance": 0.001,
        "BM_tolerance": 0.01,
        "use_vdw": False,
        "use_spin": True,
        "encut": 400,
        "k_dens": 5.0,
    },
    "sg_params": {
        "slab_thick": 4,
        "vac_thick": 20.0,
        "min_thick_A": 2.0,
        "miller": (1, 0, 0),
    },
    "sg_filter": {"method": "all"},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
    },
}

full_output = {
    "material": {"formula": "Si", "mpid": "mp-149"},
    "comp_params": {
        "functional": "PBE",
        "volume_tolerance": 0.001,
        "BM_tolerance": 0.01,
        "use_vdw": False,
        "use_spin": True,
        "encut": 400,
        "k_dens": 5.0,
        "is_metal": None,
    },
    "sg_params": {
        "slab_thick": 4,
        "vac_thick": 20.0,
        "min_thick_A": 2.0,
        "miller": (1, 0, 0),
        "max_index": None,
        "symmetrize": False,
        "primitive": True,
        "lll_reduce": True,
        "tol": 0.05,
        "max_normal_search": "max",
        "resize": True,
        "preserve_terminations": True,
        "minimize_structures": False,
        "match_ouc_lattice": True,
        "calculate_bonds": True,
        "center_slab": True,
        "max_nsites": 100,
        "filter_polar": True,
        "nn_method": "all",
        "weight": "BVS",
    },
    "sg_filter": {"method": "all", "bvs_param": None},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
        "interface_distance": "auto",
        "max_angle_tol": 0.01,
        "external_pressure": None,
    },
}


class InputCheckingTest(unittest.TestCase):
    def test_correct_test(self):
        o = check_hetero_wf_inputs(inputs_correct)
        self.assertEqual(o, full_output)

    def test_no_defaults(self):
        ic = deepcopy(inputs_correct)
        ic["material"] = {"mpid": "mp-149"}
        print(ic)
        with self.assertRaises(ValueError):
            o = check_hetero_wf_inputs(ic)

    def test_misspelled_top_level(self):
        ic = deepcopy(inputs_correct)
        ic["computational_param"] = ic.pop("comp_params")
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_misspelled_sub_level(self):
        ic = deepcopy(inputs_correct)
        ic["comp_params"]["BM_tol"] = ic["comp_params"].pop("BM_tolerance")
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_more_misspelled_sub_level(self):
        ic = deepcopy(inputs_correct)
        ic["interface_params"]["interface_gap"] = 3.4
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_correct_input_in_wrong_top_level(self):
        ic = deepcopy(inputs_correct)
        ic["sg_params"]["functional"] = "PBE"
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_no_miller_and_no_max_index(self):
        ic = deepcopy(inputs_correct)
        ic["sg_params"].pop("miller")
        with self.assertRaises(ValueError):
            o = check_hetero_wf_inputs(ic)


if __name__ == "__main__":
    unittest.main()
