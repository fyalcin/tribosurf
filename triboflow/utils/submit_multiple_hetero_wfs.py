#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 16:12:32 2021

@author: wolloch
"""

from fireworks import LaunchPad
from triboflow.workflows.main import heterogeneous_wf

materials = {
    "Al": {"mpid": "mp-134"},
    "Cu": {"mpid": "mp-30"},
    "Fe": {"mpid": "mp-13"},
}
#'C': {'mpid': 'mp-66',
#     'thick_min': 4,
#     'thick_inc': 2}}

miller_indices = ["100", "110", "111"]

default_thickness_params = {
    "thick_min": 4,
    "thick_max": 14,
    "thick_incr": 1,
}

computational_params = {
    "functional": "SCAN",
    "volume_tolerance": 0.001,
    "BM_tolerance": 0.01,
    "use_vdw": False,
    "surfene_thr": 0.01,
    "vacuum": 15,
}

interface_params = {
    "max_area": 100,
    "max_area_ratio_tol": 0.1,
    "max_length_tol": 0.025,
    "max_angle_tol": 0.01,
}


def create_triboflow_inputs(
    materials,
    miller_indices,
    default_thickness_params,
    computational_params,
    interface_params,
):
    slabs = []
    for k, v in materials.items():
        mpid = v.pop("mpid")
        for miller in miller_indices:
            m1 = {"formula": k, "mpid": mpid, "miller": miller}
            m1.update(default_thickness_params)
            m1.update(v)
            slabs.append(m1)

        inputs_list = []
        for mat_1 in slabs:
            for mat_2 in slabs:
                if not mat_1 == mat_2:
                    inputs = {
                        "material_1": mat_1,
                        "material_2": mat_2,
                        "computational_params": computational_params,
                        "interface_params": interface_params,
                    }
                    dummy_inputs = {
                        "material_1": mat_2,
                        "material_2": mat_1,
                        "computational_params": computational_params,
                        "interface_params": interface_params,
                    }
                    if dummy_inputs not in inputs_list:
                        inputs_list.append(inputs)
    return inputs_list


def submit_multiple_wfs(inputs_list):
    lpad = LaunchPad.auto_load()
    for inputs in inputs_list:
        WF = heterogeneous_wf(inputs)
        lpad.add_wf(WF)


if __name__ == "__main__":
    inputs_list = create_triboflow_inputs(
        materials=materials,
        miller_indices=miller_indices,
        computational_params=computational_params,
        interface_params=interface_params,
        default_thickness_params=default_thickness_params,
    )
    submit_multiple_wfs(inputs_list)
