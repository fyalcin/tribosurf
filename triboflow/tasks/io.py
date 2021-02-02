#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 15:26:19 2021

Useful tools to handle input/output operations on files

@author: glosi000
"""

import json

def read_json(jsonfile):
    """
    Shortcut to easily read a json file.

    Parameters
    ----------
    jsonfile : str
        Input file to be read.

    Returns
    -------
    data : dict
        Dictionary containing the data.
        
    """
    
    with open(jsonfile, 'r') as f:
        data = json.load(f)  
    return data
