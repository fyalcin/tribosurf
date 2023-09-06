import unittest
from copy import deepcopy
from triboflow.utils.check_WF_inputs import check_hetero_wf_inputs

inputs_correct = {
    "material_1": {"formula": "Cu3Sn", "mpid": "mp-149"},
    "material_2": {"formula": "SiC", "mpid": "mp-134"},
    "computational_params": {
        "functional": "SCAN",
        "volume_tolerance": 0.001,
        "BM_tolerance": 0.01,
        "use_vdw": True,
        "use_spin": True,
    },
    "sg_params_1": {
        "slab_thick": 6,
        "vac_thick": 25.0,
        "min_thick_A": 6.0,
        "miller": [1, 1, 1],
    },
    "sg_params_2": {
        "slab_thick": 6,
        "vac_thick": 25.0,
        "min_thick_A": 6.0,
        "max_index": 1,
    },
    "sg_filter_1": {"method": "bvs_min_N", "bvs_param": 2},
    "sg_filter_2": {"method": "bvs_min_N", "bvs_param": 2},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
        "external_pressure": 1,
        "max_sites": 200,
    },
    "database_params": {
        "db_file": "/home/fs71411/mwo3/FireWorks/config_VSC5_Zen3/db.json",
        "high_level": "release_test_db",
    },
}


full_output = {
    "material_1": {"formula": "Cu3Sn", "mpid": "mp-149"},
    "material_2": {"formula": "SiC", "mpid": "mp-134"},
    "computational_params": {
        "functional": "SCAN",
        "volume_tolerance": 0.001,
        "BM_tolerance": 0.01,
        "use_vdw": True,
        "use_spin": True,
        "encut": 520,
        "k_dens": 9.0,
        "is_metal": None,
    },
    "sg_params_1": {
        "slab_thick": 6,
        "vac_thick": 25.0,
        "min_thick_A": 6.0,
        "max_index": None,
        "miller": [1, 1, 1],
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
    "sg_params_2": {
        "slab_thick": 6,
        "vac_thick": 25.0,
        "min_thick_A": 6.0,
        "max_index": 1,
        "miller": None,
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
    "sg_filter_1": {"method": "bvs_min_N", "bvs_param": 2},
    "sg_filter_2": {"method": "bvs_min_N", "bvs_param": 2},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
        "external_pressure": 1,
        "max_sites": 200,
        "interface_distance": "auto",
        "max_angle_tol": 0.01,
    },
    "database_params": {
        "db_file": "/home/fs71411/mwo3/FireWorks/config_VSC5_Zen3/db.json",
        "high_level": "release_test_db",
        "bulk_coll": "bulk_data",
        "surfen_coll": "surfen_data",
        "conv_coll": "conv_data",
        "wulff_coll": "wulff_data",
    },
}


class InputCheckingTest(unittest.TestCase):
    def test_correct_test(self):
        o = check_hetero_wf_inputs(inputs_correct)
        self.assertEqual(o, full_output)

    def test_no_defaults(self):
        ic = deepcopy(inputs_correct)
        ic["material_1"] = {"mpid": "mp-149"}
        print(ic)
        with self.assertRaises(ValueError):
            o = check_hetero_wf_inputs(ic)

    def test_misspelled_top_level(self):
        ic = deepcopy(inputs_correct)
        ic["comp_params"] = ic.pop("computational_params")
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_misspelled_sub_level(self):
        ic = deepcopy(inputs_correct)
        ic["computational_params"]["BM_tol"] = ic["computational_params"].pop("BM_tolerance")
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_more_misspelled_sub_level(self):
        ic = deepcopy(inputs_correct)
        ic["interface_params"]["interface_gap"] = 3.4
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_correct_input_in_wrong_top_level(self):
        ic = deepcopy(inputs_correct)
        ic["sg_params_2"]["functional"] = "PBE"
        with self.assertRaises(KeyError):
            o = check_hetero_wf_inputs(ic)

    def test_no_miller_and_no_max_index(self):
        ic = deepcopy(inputs_correct)
        ic["sg_params_1"].pop("miller")
        with self.assertRaises(ValueError):
            o = check_hetero_wf_inputs(ic)


if __name__ == "__main__":
    unittest.main()
