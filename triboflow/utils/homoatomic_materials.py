#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 14:11:29 2021

@author: wolloch
"""
import json

import os

project_folder = os.path.dirname(__file__)


def load_homoatomic_materials(
    file_location=project_folder + "/../",
    filename="homoatomic_materials.json",
):
    with open(file_location + filename, "r") as f:
        data = json.load(f)
    return data
