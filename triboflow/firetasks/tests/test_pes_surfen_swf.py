from datetime import datetime

from fireworks import LaunchPad, Firework, Workflow

from triboflow.firetasks.run_slabs_wfs import GetCandidatesForHeteroStructure
from triboflow.firetasks.structure_manipulation import FT_MakeHeteroStructure

required_params = ["mpid_1", "mpid_2"]
optional_params = [
    "comp_params_1",
    "comp_params_2",
    "interface_params",
    "sg_params_1",
    "sg_params_2",
    "sg_filter_1",
    "sg_filter_2",
    "db_file",
    "high_level",
    "surfen_coll",
    "bulk_coll",
    "add_full_relax",
]
inputs = {
    "material_1": {"formula": "Si", "mpid": "mp-149"},
    "material_2": {"formula": "Al", "mpid": "mp-134"},
    "computational_params": {
        "functional": "PBE",
        "volume_tolerance": 0.01,
        "BM_tolerance": 0.01,
        "use_vdw": False,
        "use_spin": True,
    },
    "sg_params_1": {
        "slab_thick": 4,
        "vac_thick": 20.0,
        "min_thick_A": 3.0,
        "miller": (1, 0, 0),
    },
    "sg_params_2": {
        "slab_thick": 4,
        "vac_thick": 20.0,
        "min_thick_A": 3.0,
        "miller": (1, 0, 0),
    },
    "sg_filter_1": {"method": "bvs_min_N", "bvs_param": 1},
    "sg_filter_2": {"method": "bvs_min_N", "bvs_param": 1},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
        "external_pressure": 1,
    },
    "db_file": "/home/yalcin/config_local/db.json",
    "high_level": "surfflow_test",
}

mpid_1 = inputs["material_1"]["mpid"]
mpid_2 = inputs["material_2"]["mpid"]

wf_list = []

fw1 = Firework(
    [
        GetCandidatesForHeteroStructure(
            mpid_1=mpid_1,
            mpid_2=mpid_2,
            comp_params_1={},
            comp_params_2={},
            interface_params=inputs["interface_params"],
            sg_params_1=inputs["sg_params_1"],
            sg_params_2=inputs["sg_params_2"],
            sg_filter_1=inputs["sg_filter_1"],
            sg_filter_2=inputs["sg_filter_2"],
            db_file=inputs["db_file"],
            high_level=inputs["high_level"],
            surfen_coll="PBE.surfen_data",
            bulk_coll="PBE.bulk_data",
            add_full_relax=True,
        )
    ],
    name="Get candidates for heterostructure",
)

fw2 = Firework(
    FT_MakeHeteroStructure(
        mp_id_1=mpid_1,
        mp_id_2=mpid_2,
        interface_params=inputs["interface_params"],
        functional="PBE",
        db_file=inputs["db_file"],
        high_level=inputs["high_level"],
    ),
    name="Match the interface",
    parents=[fw1],
)

wf = Workflow([fw1, fw2], name="Test surfen workflow")
today = datetime.today().strftime("%Y-%m-%d")

lpad = LaunchPad.from_file("/home/yalcin/config_local/my_launchpad.yaml")
lpad.reset(today)
lpad.add_wf(wf)

# fworker = FWorker.from_file("/home/yalcin/config_local/my_fworker.yaml")
# rapidfire(launchpad=lpad, fworker=fworker, m_dir="/home/yalcin/scratch")
