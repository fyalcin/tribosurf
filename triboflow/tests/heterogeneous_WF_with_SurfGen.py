#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from datetime import datetime
from fireworks import LaunchPad

from triboflow.workflows.main import heterogeneous_wf

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
        "max_index": 2,
    },
    "sg_params_2": {
        "slab_thick": 4,
        "vac_thick": 20.0,
        "min_thick_A": 3.0,
        "max_index": 2,
    },
    "sg_filter_1": {"method": "bvs_min_N",
                    "bvs_param": 10},
    "sg_filter_2": {"method": "bvs_min_N",
                    "bvs_param": 10},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
        "external_pressure": 1,
    },
    "db_file": "/home/fs71411/mwo3/FireWorks/config_VSC5_Zen3/db.json",
    "high_level": "release_test_db",
}

WF = heterogeneous_wf(inputs)

lpad = LaunchPad.auto_load()
#today = datetime.today().strftime("%Y-%m-%d")
#lpad.reset(today)
lpad.add_wf(WF)

# fworker = FWorker.from_file("/home/yalcin/config_local/my_fworker.yaml")
# rapidfire(launchpad=lpad, fworker=fworker, m_dir="/home/yalcin/scratch")
