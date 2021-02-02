#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 12:35:00 2021

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__credits__ = 'This module is based on the Triboflow package, Michael Wolloch'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'January 20th, 2021'

from tribchem.workflows.REFACTORING.utils.database import NavigatorMP

def interface_name(mp_id_1, miller_1, mp_id_2, miller_2):
    """
    Return a name for an interface based on MP-IDs and miller indices.

    Parameters
    ----------
    mp_id_1 : str
        MP ID of the first material.
    miller_1 : list of int or str
        Miller indices of first material given either as list or str.
    mp_id_2 : str
        MP ID of the second material.
    miller_2 : list of int or str
        Miller indices of the second material given either as list or str.

    Returns
    -------
    name : str
        Unique name for the interface of two slabs.
    """

    nav_mp = NavigatorMP()
    f1 = nav_mp.get_property_from_mp(
        mp_id=mp_id_1, 
        properties=['pretty_formula'])['pretty_formula']
    f2 = nav_mp.get_property_from_mp(
        mp_id=mp_id_2,
        properties=['pretty_formula'])['pretty_formula']
    
    # Assign the miller index as a string
    m1 = miller_1
    m2 = miller_2
    if isinstance(miller_1, list):
        m1 = ''.join(str(s) for s in miller_1)
    if isinstance(miller_2, list):
        m2 = ''.join(str(s) for s in miller_2)

    n1 = min(f1 + m1, f2 + m2)
    n2 = max(f1 + m1, f2 + m2)
    ids = min(mp_id_1 + '_' + mp_id_2, mp_id_2 + '_' + mp_id_1)
    name = '_'.join((n1, n2, ids))
    
    return name
