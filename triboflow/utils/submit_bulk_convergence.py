#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 16:12:32 2021

@author: wolloch
"""

from fireworks import LaunchPad
from triboflow.workflows.main import optimize_bulk_wf
from triboflow.utils.homoatomic_materials import load_homoatomic_materials

computational_params = {
    "functional": "SCAN",
    "volume_tolerance": 0.001,
    "BM_tolerance": 0.01,
    "use_vdw": True,
}
db_file = "/home/fs71411/mwo3/FireWorks/config_VSC5_Zen3/db.json"
high_level = "release_test_db"


def submit_multiple_wfs(workflow_list):
    lpad = LaunchPad.auto_load()
    for wf in workflow_list:
        lpad.add_wf(wf)


if __name__ == "__main__":
    materials_dict = load_homoatomic_materials()

    # materials_list = []
    # for formula in ['Al', 'C', 'Si', 'Ge', 'Cu', 'Ag',
    #                 'Au', 'Ni', 'Fe', 'Ti', 'Co']:
    #     mpid = materials_dict[formula]['mpids'][materials_dict[formula]['default']]
    #     materials_list.append({'formula': formula, 'mpid': mpid})

    materials_list = [
        {"formula": "Al2O3", "mpid": "mp-1143"},
        {"formula": "Cu3Sn", "mpid": "mp-13138"},
    ]

    workflow_list = []
    for mat in materials_list:
        inputs = {
            "material": mat,
            "computational_params": computational_params,
            "db_file": db_file,
            "high_level": high_level,
        }
        workflow_list.append(optimize_bulk_wf(inputs))
    submit_multiple_wfs(workflow_list)
