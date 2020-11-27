#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 09:50:38 2020

@author: mwo
"""

import numpy as np
import math as m
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from pymatgen.core.surface import Slab
from pymatgen.core.operations import SymmOp
from triboflow.phys.high_symmetry import GetSlabHS, GetInterfaceHS, \
    PBC_HSPoints, RemoveDuplicatesFromHSDicts
from triboflow.phys.potential_energy_surface import UnfoldPES
from triboflow.utils.database import GetDBJSON, GetInterfaceFromDB
from triboflow.utils.plot_tools import Plot_SlabHS

def test_function(x,y):
    r = (x**2 + y**2)**0.5
    z = np.cos(2*r) * np.exp(-0.5*abs(r))
    return z

db_file = GetDBJSON()
functional = 'PBE'
name = 'FeRh001_MgO001_mp-1265_mp-1918'

# =============================================================================
# rand_points = np.random.rand(500,2)*4*m.pi - 2*m.pi
# x = rand_points[:,0]
# y = rand_points[:,1]
# fig = plt.figure(figsize=(7, 7), dpi=100)
# ax = fig.add_subplot(111)
# ax.set_aspect('equal')
# plt.plot(x,y, 'ro' )
# z = test_function(x, y)
# 
# # fig = plt.figure()
# # ax = Axes3D(fig)
# 
# # ax.scatter(x, y, z)
# 
# rbf = Rbf(x, y, z)
# 
# grid_x = np.arange(min(x), max(x), 0.05)
# grid_y = np.arange(min(y), max(y), 0.05)
# 
# X, Y = np.meshgrid(grid_x, grid_y)
# Z = rbf(X, Y)
# 
# #plt.contourf(X, Y, Z)
# 
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# 
# fig.colorbar(surf, shrink=0.5, aspect=5)
# =============================================================================


interface_dict = GetInterfaceFromDB(name, db_file, functional)
c_u = interface_dict['PES']['high_symmetry_points']['combined_unique']
c_a = interface_dict['PES']['high_symmetry_points']['combined_all']
E_l = interface_dict['PES']['high_symmetry_points']['energy_list']

top_aligned = Slab.from_dict(interface_dict['top_aligned'])
bottom_aligned = Slab.from_dict(interface_dict['bottom_aligned'])
        
#top slab needs to be mirrored to find the high symmetry points at the
#interface.
mirror = SymmOp.reflection(normal=[0,0,1], origin=[0, 0, 0])
flipped_top = top_aligned.copy()
flipped_top.apply_operation(mirror)
top_hsp_unique, top_hsp_all = GetSlabHS(flipped_top)
        
bottom_hsp_unique, bottom_hsp_all = GetSlabHS(bottom_aligned)
        
cell = bottom_aligned.lattice.matrix
        
hsp_unique = GetInterfaceHS(bottom_hsp_unique, top_hsp_unique, cell)
hsp_all = GetInterfaceHS(bottom_hsp_all, top_hsp_all, cell)
        
c_hsp_u, c_hsp_a = RemoveDuplicatesFromHSDicts(hsp_unique,
                                               hsp_all,
                                               decimals=5)

Plot_SlabHS(top_hsp_unique, top_aligned)
Plot_SlabHS(top_hsp_all, top_aligned)
Plot_SlabHS(bottom_hsp_unique, bottom_aligned)
Plot_SlabHS(bottom_hsp_all, bottom_aligned)
Plot_SlabHS(c_hsp_u, top_aligned)
Plot_SlabHS(c_hsp_a, top_aligned)


double_keys = {}
for key, value in c_hsp_a.items():
    for point in value:
        double_keys[key] = []
        for key_2, value_2 in c_hsp_a.copy().items():
            if point in value_2 and key != key_2:
                double_keys[key].append(key_2)
                
E_reduced = []
for l in E_l:
    if l[0] in c_hsp_u.keys():
        E_reduced.append(l)
        
c_hsp_all_reduced = {}
for key, value in c_hsp_a.items():
     array = np.unique(np.array(value), axis=0)
     c_hsp_all_reduced[key] = array.tolist()

E_list, data = UnfoldPES(c_hsp_all_reduced, E_reduced)