#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Test PES/shear strength modules, using test points
Created on Mon Oct  26 11:20:07 2020

@author: gl
"""

import numpy as np
from PES_functions import *
from utility_functions import *

pes_array = np.genfromtxt('pes.dat')
pes_array = np.unique(pes_array, axis=0)

StaticTribo(interface, hs, E)