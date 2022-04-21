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
import matplotlib.pyplot as plt

from PIL import Image
from io import BytesIO
from scipy.interpolate import RBFInterpolator

class PESGenerator():
    def __init__(self,
                 points_per_angstrom = 50,
                 interpolation_kernel = 'linear',
                 plot_hs_points = False,
                 plot_unit_cell = True,
                 plotting_ratio = 1.0,
                 normalize_minimum = True,
                 nr_of_contours = 30,
                 fig_title = 'PES',
                 fig_type = 'png',
                 plot_path = None):
        """
        Class for generating PES interpolation and plots for interfaces.
        
        Given an interface, and a couple of dictionaries, see the __call__
        method, a smooth PES is constructed using radial basis function
        interpolation. The interpolation object, as well as a figure of the
        PES are provided as properties of the class once the class has been
        called as a funcion on the interface object and the necessary dicts
        with high-symmetry point and energy data.
        
        
        Important properties of the class:
            .rbf contains the scipy.interpolate._rbfinterp.RBFInterpolator
            .PES_fig contains the matplotlib.figure.Figure of the PES
            .PES_as_bytes contains the PES as a bytes object.
        
        
        
        Parameters
        ----------
        points_per_angstrom : int, optional
            Points in the interpolation meshgrid per Angstrom.
            The default is 50.
        interpolation_kernel : str, optional
            Option for the RBFInterpolator. Possibilities are: 
                - 'linear'               : ``-r``
                - 'thin_plate_spline'    : ``r**2 * log(r)``
                - 'cubic'                : ``r**3``
                - 'quintic'              : ``-r**5``
                - 'multiquadric'         : ``-sqrt(1 + r**2)``
                - 'inverse_multiquadric' : ``1/sqrt(1 + r**2)``
                - 'inverse_quadratic'    : ``1/(1 + r**2)``
                - 'gaussian'             : ``exp(-r**2)``
            The default is 'linear' and should probably not be changed to
            keep proper symmetry!
        plot_hs_points : bool or str, optional
            Plot 'all' or 'unique' high symmetry point if those strings are
            passed. If not, no points are plotted. The default is False.
        plot_unit_cell : bool, optional
            Plot the unit cell in the PES. The default is True.
        plotting_ratio : float or None, optional
            Aspect ratio of the PES figure. Can be set to None in which case
            the smallest rectangle fitting the unit cell is plotted.
            The default is 1.0.
        normalize_minimum : bool, optional
            Set the minimum of the interpolation to 0 and make all other
            energies positive. The default is True.
        nr_of_contours : int, optional
            how many contours there are in the PES figure. The default is 30.
        fig_title : str, optional
            Title of the PES figure. The default is 'PES'.
        fig_type : str, optional
            Selects the format of the plot, e.g. 'pdf', or 'svg'.
            The default is 'png'.
        plot_path : str or None, optional
            If a string is set, the PES will be saved at this path. USE A
            TRAILING SLASH! E.g.: /home/user/test/PES/
            The default is None.

        """
        self.ppA = points_per_angstrom
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
        """
        Interpolate a PES for a given interface and prepare a plot of it.

        Parameters
        ----------
        interface : pymatgen.core.interface.Interface
            A pymatgen interface object
        energies_dict : dict
            Adhesion energies for each group of high symmetry points.
        all_shifts_dict : dict
            All lateral shifts grouped by high-symmetry group.
        unique_shifts_dict : dict
            Unique lateral shifts grouped by high-symmetry group.
        group_names_dict : dict, optional
            Assigning meaningful shift names to the group names.
            The default is None.

        """
        self.interface = interface.copy()
        self.energies_dict = energies_dict.copy()
        self.all_shifts = all_shifts_dict.copy()
        self.unique_shifts = unique_shifts_dict.copy()
        if group_names_dict:
            self.group_names_dict = group_names_dict.copy()
        else:
            self.group_names_dict = None
        
        self.__get_limits_and_multiples()

        X, Y, Z = self.__interpolate_on_grid()
        
        self.__plot_grid(X, Y, Z)
        self.PES_as_bytes = self.__get_pes_as_bytes(self.PES_fig)
    
        
    def __get_limits_and_multiples(self):
        """
        Compute the expansion of the input data and plot limits using unit cell
        and aspect ratio.

        """
        self.__get_plotting_rectangle()
        self.__set_meshgrid_limits()
        self.__get_mult_limits()
        
    def __set_meshgrid_limits(self):
        """
        Set plot and meshgrid limits according to aspect ratio

        """
        if self.plotting_ratio:
            if self.width/self.height < self.plotting_ratio:
                xlim = self.height*self.plotting_ratio
                ylim = self.height
            else:
                xlim = self.width
                ylim = self.width/self.plotting_ratio
        else:
            xlim = self.width
            ylim = self.height
        self.xlim = xlim
        self.ylim = ylim
        
    def __get_mult_limits(self):
        """
        Set the limits (mostly generously) of data replication to get a proper
        interpolation and avoid edge effects.
        """
        a = self.interface.lattice.matrix[0,:2]
        b = self.interface.lattice.matrix[1,:2]
        ab = a+b
        max_a = max(a[0], ab[0])
        max_b = max(b[1], ab[1])
        self.xmult = np.ceil(self.xlim/max_a)+1
        self.ymult = np.ceil(self.ylim/max_b)+1
        
    def __make_energies_list(self):
        """
        Create a list of points with their respective PES energies.
        
        [[x1, y1, E1],
         [x2, y2, E2],
         :
         [xn, yn, Em]]

        """
        energy_list = []
        for group, coordinates in self.all_shifts.items():
            for coord in coordinates:
                energy_list.append(coord + [self.energies_dict[group]])
        self.unit_cell_energies = energy_list
        
    def __extend_energies_list(self):
        """
        Extend the input data from the unit cell to unit cells around it to
        avoid edge effects in interpolation and allow to conform to chosen
        aspect ratios.

        """
        extended_energy_list = []
        xrange = np.arange(-self.xmult, self.xmult, 1.0)
        yrange = np.arange(-1, self.ymult, 1.0)
        for x in xrange:
            for y in yrange:
                for entry in self.unit_cell_energies:
                    extended_energy_list.append([entry[0]+x,
                                                 entry[1]+y,
                                                 entry[2]])
        self.extended_energies = self.__from_frac_to_cart(
            np.asarray(extended_energy_list))
        
    def __get_rbf(self):
        """
        Create a RBF function from the replicated input data.
        """
        self.__make_energies_list()
        self.__extend_energies_list()
        self.rbf = RBFInterpolator(self.extended_energies[:, :2],
                        self.extended_energies[:, 2],
                        kernel=self.interpolation_kernel)
        
    def __from_frac_to_cart(self, array):
        """
        Transform fractional to cartesian coordinates

        Parameters
        ----------
        array : list or np.ndarray
            Array of points in fractional coordinates

        Returns
        -------
        numpy.ndarray
            Input array transformed to cartesian coordinates

        """
        m = self.interface.lattice.matrix[:2,:2]
        if np.asarray(array).ndim == 1:
            array = [array]
        if len(array[0]) == 3:
            m = np.vstack((np.hstack((m,np.zeros((2,1)))),np.asarray([0,0,1])))
        return np.dot(array,m)
        
                      
            
    def __get_plotting_rectangle(self):
        """
        Get a rectangle that is just large enough to contain the unit cell.

        """
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
        """
        Plot the PES in a matplotlib figure adding unit cell and high symmetry
        points if wanted.

        Parameters
        ----------
        X : nump.yndarray
            X part of meshgrid
        Y : nump.yndarray
            Y part of meshgrid
        Z : nump.yndarray
            Interpolation ready for plotting
        """
        levels = np.linspace(np.amin(Z), np.amax(Z), self.contours)
        fig, ax = self.__get_fig_and_ax()
        zt1 = plt.contourf(X, Y, Z, levels, cmap=plt.cm.RdYlBu_r)
        ticks = self.__colorbar_ticks(Z)

        cbar = plt.colorbar(zt1, ax=ax, orientation='vertical',
                             ticks=ticks,
                             format='%.3f', shrink=0.7)
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label(r'$\rm{E_{adh}} [\rm J/ \rm m^2]$', rotation=270, labelpad=25,
                        fontsize=18, family='sans-serif')
        if self.plot_unit_cell:
            self.__plot_unit_cell(ax)
        
        plt.xlabel(r"x [$\rm\AA$]", fontsize=20, family='sans-serif')
        plt.ylabel(r"y [$\rm\AA$]", fontsize=20, family='sans-serif')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim = ((0,self.xlim))
        plt.ylim = ((0,self.ylim))
        ax.set_title(self.fig_name, fontsize = 24, family='sans-serif', pad=10)
        if self.plot_hs_points == 'unique':
            self.__plot_hs_points(self.unique_shifts, fig, self.group_names_dict)
        elif self.plot_hs_points == 'all':
            self.__plot_hs_points(self.all_shifts, fig, self.group_names_dict)
        
        
        fig.set_tight_layout(True)
        plt.tight_layout()
        self.PES_fig = fig
        
        if self.plot_path:
            plt.savefig(self.plot_path+self.fig_name+'.'+self.fig_type,
                        dpi=300, bbox_inches='tight')
              
    def __evaluate_on_grid(self, X, Y):
        """
        Evalueate the RBF on a meshgrid

        Parameters
        ----------
        X : nump.yndarray
            X part of meshgrid
        Y : nump.yndarray
            Y part of meshgrid

        Returns
        -------
        Z : nump.yndarray
            Interpolation ready for plotting.
        """
        xy = np.vstack([X.ravel(), Y.ravel()]).T
        z = self.rbf(xy)
        Z = np.reshape(z, X.shape, order='C')
        if self.normalize:
            Z = Z - min(Z.ravel())
        return Z
    
    def __get_grid(self, xmax, ymax):
        """
        Return a meshgrid for a (0,xmax*xmult) (0,ymax*ymult) rectangle

        Parameters
        ----------
        xmax : float
            width of the rectangle for interpolation
        ymax : float
            height of the rectangle for interpolation

        Returns
        -------
        X : nump.yndarray
            X part of meshgrid
        Y : nump.yndarray
            Y part of meshgrid
        """
        nr_pts_x = int(xmax*self.ppA)
        nr_pts_y = int(ymax*self.ppA)
        grid_x = np.linspace(0, xmax, nr_pts_x)
        grid_y = np.linspace(0, ymax, nr_pts_y)
        X, Y = np.meshgrid(grid_x, grid_y)
        return X, Y
        
    def __colorbar_ticks(self, Z, nr_of_ticks=10):
        """
        Define tick marks for the colorbar. Make sure the lowest and highest
        values are included

        Parameters
        ----------
        Z : numpy.ndarray
            Interpolated PES values
        nr_of_ticks : int, optional
            Number of tick marks. The default is 8.

        Returns
        -------
        numpy.ndarray
            Array of tick marks for the colorbar

        """
        low = min(Z.ravel())
        high = max(Z.ravel())
        return np.linspace(low, high, nr_of_ticks).tolist()
    
    def __plot_unit_cell(self, ax):
        """
        Adds the unit cell to the plot

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            Subplot axes
        """
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
        """
        Generate a matplotlib figure and corresponding axes.
        
        Figure size depends on the aspect ratio and if high symmetry points
        should be plotted.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure
        ax : matplotlib.axes._subplots.AxesSubplot
            Subplot axes

        """
        add_height = self.ylim*0.10 if self.plot_hs_points else 0.0
        if self.plotting_ratio:
            fig = plt.figure(figsize=(self.xlim,
                                      self.xlim/self.plotting_ratio+add_height),
                             dpi=300)
        else:
            fig = plt.figure(figsize=(self.xlim, self.ylim+add_height),
                             dpi=300)
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        return fig, ax
            
    def __plot_hs_points(self, hs_points_dict, fig, group_names_dict=None):
        """
        Plot high symmetry points including a legend on the PES

        Parameters
        ----------
        hs_points_dict : dict
            Dictionary holding group names and 2D point arrays
        fig : matplotlib.figure.Figure
            PES figure
        group_names_dict : dict, optional
            Dictionary mapping the group names to meaningful labels, e.g.:
            'ontop_1-hollow_2'. The default is None.

        Returns
        -------
        None.

        """
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
                       columnspacing=1.0,
                       labelspacing=0.3,
                       #ncol=math.ceil(len(hs_points_dict)/2),
                       ncol=3,
                       fontsize=15)
    
    def __get_pes_as_bytes(self,fig):
        """
        Transform a figure into a bytes object to store in the MongoDB database.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure to be converted

        Returns
        -------
        bytes
            Figure in bytes

        """
        im = Image.frombytes('RGB', 
                             fig.canvas.get_width_height(),
                             fig.canvas.tostring_rgb())
        image_bytes = BytesIO()
        im.save(image_bytes, format='png')
        return image_bytes.getvalue()
    
    def __interpolate_on_grid(self):
        """
        Generate RBF interpolation, get a meshgrid and evaluate the RBF there.
        
        Interpolation is done initially on a rectangular meshgrid that just
        fits the unit cell. Edge effects are circumvented by replicating the
        input date laterally.

        Returns
        -------
        X : numpy.ndarray
            X part of the meshgrid
        Y : numpy.ndarray
            Y part of the meshgrid
        Z : numpy.ndarray
            Interpolation on the meshgrid ready for plotting or replicating

        """
        self.__get_rbf()
        
        X, Y = self.__get_grid(self.xlim, self.ylim)
        Z = self.__evaluate_on_grid(X, Y)
        return X, Y, Z
    