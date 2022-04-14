#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 14 2022

Python class to interpolate and plot the Potential Energy Surface (PES) of an interface.

"""

__author__ = 'Michael Wolloch'
__copyright__ = 'Copyright 2022, Michael Wolloch, HIT, University of Vienna'
__credits__ = 'Code partly inspired by a previous version of M. Wolloch and Gabriele Losi'
__contact__ = 'michael.wolloch@univie.ac.at'
__date__ = 'April 14th, 2022'


import numpy as np
import math
import matplotlib.pyplot as plt

from PIL import Image
from io import BytesIO
from scipy.interpolate import RBFInterpolator

class PESGenerator():
    def __init__(self,
                 point_density = 50,
                 interpolation_kernel = 'linear',
                 plot_hs_points = 'all',
                 plot_unit_cell = True,
                 plotting_ratio = 16/9,
                 normalize_minimum = False,
                 nr_of_contours = 30,
                 fig_title = 'PES',
                 fig_type = 'png',
                 plot_path = None):
        self.point_density = point_density
        self.interpolation_kernel = interpolation_kernel
        self.plotting_ratio = plotting_ratio
        self.normalize = normalize_minimum
        self.contours = nr_of_contours
        self.plot_unit_cell = plot_unit_cell
        self.plot_hs_points = plot_hs_points
        self.fig_name = fig_title
        self.fig_type = fig_type
        self.plot_path = plot_path
        
    def __call__(self,
                 interface,
                 energies_dict,
                 all_shifts_dict,
                 unique_shifts_dict,
                 group_names_dict=None):
        
        self.interface = interface.copy()
        self.energies_dict = energies_dict.copy()
        self.all_shifts = all_shifts_dict.copy()
        self.unique_shifts = unique_shifts_dict.copy()
        if group_names_dict:
            self.group_names_dict = group_names_dict.copy()
        
        X, Y, Z = self.__interpolate_on_grid()
        X_ext, Y_ext, Z_ext = self.__extend_plotting_range(Z)
        
        self.__plot_grid(X_ext, Y_ext, Z_ext)
        self.plot_in_bytes = self.__get_pes_as_bytes(self.PES_fig)
        
        
    def __make_energies_list(self):
        
        energy_list = []
        for group, coordinates in self.all_shifts.items():
            for coord in coordinates:
                energy_list.append(coord + [self.energies_dict[group]])
        self.unit_cell_energies = energy_list
        
    def __extend_energies_list(self):
        extended_energy_list = []
        n_extend = 1.0
        xrange = yrange = np.arange(-n_extend, n_extend+1.0, 1.0)
        for x in xrange:
            for y in yrange:
                for entry in self.unit_cell_energies:
                    extended_energy_list.append([entry[0]+x,
                                                 entry[1]+y,
                                                 entry[2]])
        self.extended_energies = self.__from_frac_to_cart(
            np.asarray(extended_energy_list))
        
    def __get_rbf(self):
        self.__make_energies_list()
        self.__extend_energies_list()
        self.rbf = RBFInterpolator(self.extended_energies[:, :2],
                        self.extended_energies[:, 2],
                        kernel=self.interpolation_kernel)
        
    def __from_frac_to_cart(self, array):
        m = self.interface.lattice.matrix[:2,:2]
        if np.asarray(array).ndim == 1:
            array = [array]
        if len(array[0]) == 3:
            m = np.vstack((np.hstack((m,np.zeros((2,1)))),np.asarray([0,0,1])))
        return np.dot(array,m)
            
    def __get_plotting_rectangle(self):
        a = self.interface.lattice.matrix[0,:2]
        b = self.interface.lattice.matrix[1,:2]
        ab = a+b
        min_x = min(0, a[0], b[0], ab[0])
        max_x = max(a[0], b[0], ab[0])
        min_y = min(0, a[1], b[1], ab[1])
        max_y = max(a[1], b[1], ab[1])
        width = max_x - min_x
        height = max_y - min_y
        self.shift_x = min_x
        self.shift_y = min_y
        self.width = width
        self.height = height
        
    def __plot_grid(self, X, Y, Z):
        levels = np.linspace(np.amin(Z), np.amax(Z), self.contours)
        fig, ax = self.__get_fig_and_ax()
        zt1 = plt.contourf(X, Y, Z, levels, cmap=plt.cm.RdYlBu_r)
        ticks = self.__colorbar_ticks(Z)

        cbar = plt.colorbar(zt1, ax=ax, orientation='vertical',
                             ticks=ticks, format='%.3f', shrink=0.7)
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label(r'$\rm{E_{adh}} [\rm J/ \rm m^2]$', rotation=270, labelpad=25,
                        fontsize=18, family='sans-serif')
        if self.plot_unit_cell:
            self.__plot_unit_cell(ax)
        
        plt.xlabel(r"x [$\rm\AA$]", fontsize=20, family='sans-serif')
        plt.ylabel(r"y [$\rm\AA$]", fontsize=20, family='sans-serif')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        self.__set_plot_limits(plt)
        if self.plot_hs_points == 'unique':
            self.__plot_hs_points(self.unique_shifts, plt, fig, self.group_names_dict)
        elif self.plot_hs_points == 'all':
            self.__plot_hs_points(self.all_shifts, plt, fig, self.group_names_dict)
        
        ax.set_title(self.fig_name, fontsize = 24, family='sans-serif')
        fig.set_tight_layout(True)
        plt.tight_layout()
        self.PES_fig = fig
        
        if self.plot_path:
            plt.savefig(self.plot_path+self.fig_name+'.'+self.fig_type,
                        dpi=500, bbox_inches='tight')
        

        
    def __evaluate_on_grid(self, X, Y):
        xy = np.vstack([X.ravel(), Y.ravel()]).T
        z = self.rbf(xy)
        Z = np.reshape(z, X.shape, order='C')
        if self.normalize:
            Z = Z - min(Z.ravel())
        return Z
    
    def __get_grid(self, xmax, ymax, xmult=1, ymult=1):
        nr_pts_x = int(xmax*self.point_density)*xmult
        nr_pts_y = int(ymax*self.point_density)*ymult
        grid_x = np.linspace(0, xmax*xmult, nr_pts_x)
        grid_y = np.linspace(0, ymax*ymult, nr_pts_y)
        X, Y = np.meshgrid(grid_x, grid_y)
        return X, Y
            
    def __find_rect_multiplicator(self, width, height):
        f = self.plotting_ratio / (width/height)
        if f > 1:
            return int(np.ceil(f)), 1
        elif f < 1:
            return 1, int(np.ceil(1/f))
        else:
            return 1, 1
        
    def __colorbar_ticks(self, Z, nr_of_ticks=8):
        low = min(Z.ravel())
        high = max(Z.ravel())
        return np.linspace(low, high, nr_of_ticks).tolist()
    
    def __plot_unit_cell(self, ax):
        a = self.interface.lattice.matrix[0]
        b = self.interface.lattice.matrix[1]
        import matplotlib.patches as patches
        x = [0, a[0], a[0] + b[0], b[0]]
        x_shifted = [i-self.shift_x for i in x]
        y = [0, 0, b[1], a[1] + b[1], b[1]]
        y_shifted = [i-self.shift_y for i in y]
        ax.add_patch(patches.Polygon(xy=list(zip(x_shifted, y_shifted)),
                                     fill=False, lw=2))
    
    def __get_fig_and_ax(self):
        
        add_height = self.height*0.25 if self.plot_hs_points else 0.0
        if self.plotting_ratio:
            fig = plt.figure(figsize=(self.width,
                                      self.width/self.plotting_ratio+add_height),
                             dpi=600)
        else:
            fig = plt.figure(figsize=(self.width, self.height+add_height),
                             dpi=600)
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        return fig, ax
    
    def __set_plot_limits(self, plt):
        if self.plotting_ratio:
            x_mult, y_mult = self.__get_x_and_y_mult()
            if x_mult/y_mult < self.plotting_ratio:
                plt.xlim((0, self.height*y_mult*self.plotting_ratio))
                plt.ylim((0, self.height*y_mult))
            else:
                plt.xlim((0, self.width*x_mult))
                plt.ylim((0, self.width*x_mult/self.plotting_ratio))
        else:
            plt.xlim((0, self.width))
            plt.ylim((0, self.height))
            
    def __plot_hs_points(self, hs_points_dict, plt, fig, group_names_dict=None):
        for key, shifts in hs_points_dict.items():
            if self.group_names_dict:
                label = group_names_dict[key]
            else:
                label = key
            shifts = self.__from_frac_to_cart(shifts)
            plt.plot(shifts[:,0]-self.shift_x,
                     shifts[:,1]-self.shift_y,
                     label=label,
                     marker = 'o',
                     linestyle = '',
                     markeredgecolor='k',
                     markersize=8,
                     markeredgewidth=0.5,
                     zorder=1000.0)
            fig.legend(loc='upper left',
                       handletextpad=0.01,
                       #ncol=math.ceil(len(hs_points_dict)/2),
                       ncol=3,
                       fontsize=15)
    
    def __get_pes_as_bytes(self,fig):
        im = Image.frombytes('RGB', 
                             fig.canvas.get_width_height(),
                             fig.canvas.tostring_rgb())
        image_bytes = BytesIO()
        im.save(image_bytes, format='png')
        return image_bytes.getvalue()
    
    def __get_x_and_y_mult(self):
        if self.plotting_ratio:
            x_mult, y_mult = self.__find_rect_multiplicator(self.width,
                                                           self.height)
        else:
            x_mult = y_mult = 1
        return x_mult, y_mult
    
    def __interpolate_on_grid(self):
        self.__get_rbf()
        self.__get_plotting_rectangle()
        
        X, Y = self.__get_grid(self.width, self.height)
        Z = self.__evaluate_on_grid(X, Y)
        return X, Y, Z
    
    def __extend_plotting_range(self, Z):
        x_mult, y_mult = self.__get_x_and_y_mult()
        X_ext, Y_ext = self.__get_grid(self.width, self.height, x_mult+1, y_mult+1)
        Z_ext = np.block([[Z for i in range(x_mult+1)] for j in range(y_mult+1)])
        return X_ext, Y_ext, Z_ext
    