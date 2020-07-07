#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:10:30 2020

@author: mwo
"""

from triboflow.helper_functions import GetDB


db = GetDB()

db.drop_collection('locpot_fs.chunks')
db.drop_collection('locpot_fs.files')
db.drop_collection('aeccar2_fs.chunks')
db.drop_collection('aeccar2_fs.files')
db.drop_collection('aeccar0_fs.chunks')
db.drop_collection('aeccar0_fs.files')
db.drop_collection('chgcar_fs.chunks')
db.drop_collection('chgcar_fs.files')
