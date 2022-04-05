#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 18:12:33 2020

Python functions to get the Potential Energy Surface (PES) of an interface.

The module contains the following functions:

    - get_pes
    - remove_duplicates
    - unfold_pes
    - plot_pes

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna
    Code readapted from our past homogeneous workflow, MIT license,
    https://github.com/mcrighi/interface-workflow

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__credits__ = 'Code readapted from our past homogeneous workflow, MIT license, https://github.com/mcrighi/interface-workflow,'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 8th, 2021'

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import Rbf, RBFInterpolator

from triboflow.utils.phys_tools import replicate_points, generate_uniform_grid,\
    orthorombize, pbc_coordinates

# =============================================================================
# EVALUATION OF THE PES - MAIN
# =============================================================================


class PESGenerator():
    def __init__(self,
                 point_density=50,
                 interpolation_kernel='cubic',
                 plot_unique_points=True,
                 plot_unit_cell=True,
                 plotting_ratio=None,
                 normalize_minimum = True,
                 nr_of_contours = 50):
        self.point_density = point_density
        self.interpolation_kernel = interpolation_kernel
        self.plotting_ratio = plotting_ratio
        self.normalize = normalize_minimum
        self.contours = nr_of_contours
        self.plot_unit_cell = plot_unit_cell
        
    def __call__(self,
                 interface,
                 energies_dict,
                 all_shifts_dict,
                 unique_shifts_dict):
        
        self.interface = interface.copy()
        self.energies_dict = energies_dict.copy()
        self.all_shifts = all_shifts_dict.copy()
        self.unique_shifts = unique_shifts_dict.copy()
        if self.normalize:
            self.energy_offset = min(self.energies_dict.values())
        else:
            self.energy_offset = 0.0
        self._get_rbf()
        self._get_plotting_rectangle()
        X, Y = self._get_grid(self.width, self.height)
        Z = self._evaluate_on_grid(X, Y)
        self._plot_grid(X, Y, Z)
        
    def __make_energies_list(self):
        
        energy_list = []
        for group, coordinates in self.all_shifts.items():
            for coord in coordinates:
                energy_list.append(coord + [self.energies_dict[group]-
                                            self.energy_offset])
        self.unit_cell_energies = energy_list
        
    def __extend_energies_list(self):
        extended_energy_list = []
        xrange = yrange = [-1.0, 0.0, 1.0]
        for x in xrange:
            for y in yrange:
                for entry in self.unit_cell_energies:
                    extended_energy_list.append([entry[0]+x,
                                                 entry[1]+y,
                                                 entry[2]])
        self.extended_energies = self._from_frac_to_cart(
            np.asarray(extended_energy_list))
        
    def _get_rbf(self):
        self.__make_energies_list()
        self.__extend_energies_list()
        # self.rbf = RBFInterpolator(self.extended_energies[:, :2],
        #                self.extended_energies[:, 2],
        #                kernel=self.interpolation_kernel)
        self.rbf = Rbf(self.extended_energies[:, 0],
                       self.extended_energies[:, 1],
                       self.extended_energies[:, 2],
                       function=self.interpolation_kernel)
        
    def _from_frac_to_cart(self, array):
        m = self.interface.lattice.matrix[:2,:2]
        if len(array[0]) == 3:
            m = np.vstack((np.hstack((m,np.zeros((2,1)))),np.asarray([0,0,1])))
        return np.dot(array,m)
            
    def _get_plotting_rectangle(self):
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
        if self.plotting_ratio:
            x_mult, y_mult = self._find_rect_multiplicator(width, height)
            width = x_mult*width
            height = y_mult*height
        self.width = width
        self.height = height
        
    def _plot_grid(self, X, Y, Z):
        levels = np.linspace(np.amin(Z), np.amax(Z), self.contours)
        fig = plt.figure(figsize=(self.width, self.height), dpi=600)
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        zt1 = plt.contourf(X, Y, Z, levels, cmap=plt.cm.RdYlBu_r)
        cbar1 = plt.colorbar(zt1, ax=ax, orientation='vertical')
        cbar1.set_label(r'$E_{adh} \left(J/m^2\right)$', rotation=270, labelpad=20,
                        fontsize=15, family='sans-serif')
        if self.plot_unit_cell:
            a = self.interface.lattice.matrix[0]
            b = self.interface.lattice.matrix[1]
            import matplotlib.patches as patches
            x = [0, a[0], a[0] + b[0], b[0]]
            x_shifted = [i-self.shift_x for i in x]
            y = [0, 0, b[1], a[1] + b[1], b[1]]
            y_shifted = [i-self.shift_y for i in y]
            ax.add_patch(patches.Polygon(xy=list(zip(x_shifted, y_shifted)), fill=False, lw=2))
        plt.savefig('PES_test.svg', dpi=600)
        
    def _evaluate_on_grid(self, X, Y):
        
        Z = self.rbf(X, Y)
        return Z
    
    def _get_grid(self, xmax, ymax):
        grid_x = np.arange(0, xmax, 1.0/self.point_density)
        grid_y = np.arange(0, ymax, 1.0/self.point_density)
        X, Y = np.meshgrid(grid_x, grid_y)
        return X, Y
            
    def _find_rect_multiplicator(self, width, height):
        f = self.plotting_ratio / (width/height)
        if f > 1:
            return 1.0, np.ceil(f)
        elif f < 1:
            return np.ceil(1/f), 1.0
        else:
            return 1.0, 1.0
        

def get_pes(hs_all, E, cell, to_fig=None, point_density=20):
    """Interpolate the PES using a list of high symmetry points and energies.
    
    Main function to get the Potential Energy Surface (PES) for an interface. 
    The points are replicated to span a 3x3 lattice cell and are interpolated
    by using Radial Basis Functions (cubic function).
    In the output data the energy is normalized so that the absolute minimum
    is 0. Furthermore it is made sure that the lateral point are inside the
    unit cell before they are replicated.
    
    Parameters
    ----------        
    hs : dict
        Unfolded HS points of the interface, covering all the surface.
        
    E : dict
        Contains the energy calculated for each unique surfacial HS site.
    
    cell : numpy.ndarray
        Vectors of the lattice cell of the interface.
        
    #to_fig : string, optional CURRENTLY NOT IMPLEMENTED
    #    Name of the image that you want to save, it will be: 'name_fig'+'.pdf' 
    #    Suggested name: name = 'PES_' + 'Name of the interface'.
    #    The default is None and no image is saved.

    Returns
    -------
    rbf : scipy.interpolate.rbf.Rbf
        Object containing the information of the interpolation of the potential
        energy. Call it on a set of [x, y] coordinates to obtain the energy 
        values for those points. Usage: rbf([x, y])
        
    pes_dict : dict
        Dictionary containing the points and the associated energies divided
        by type. Useful to keep track of the type of the HS sites and their
        energies. To be saved to DB (?)
    
    pes_data : numpy.ndarray
        The entire set of HS points covering the interface with the 
        corresponding energies.
        Format:
            x[0]  y[0]  E[0]
            x[1]  y[1]  E[1]
             .     .     .
             .     .     .
             .     .     .
             
    """

    # Unfold the PES points
    E_list, data = unfold_pes(hs_all, E)
    
    #making sure points are not represented twice by ensuring rows in data are unique
    data = remove_duplicates(data)
    #print(len(data))
    
    #make sure that the x and y coordinates are inside the unit cell.
    # x_y_insideCell = pbc_coordinates(data[:, :2],
    #                                  cell,
    #                                  to_array=True)
    # data[:, :2] = x_y_insideCell
        
    # Normalize the minimum to 0
    data[:,2] = (data[:,2]-min(data[:,2]))
    
    # Interpolate the data with Radial Basis Function
    data_rep = replicate_points(data, cell, replicate_of=(3, 3) )
    rbf = Rbf(data_rep[:, 0], data_rep[:, 1], data_rep[:, 2], function='cubic')
    # rbf = RBFInterpolator(data_rep[:,:2], data_rep[:, 2], kernel='cubic')

    xrange, yrange = sum(abs(cell[:2, 0])), sum(abs(cell[:2, 1]))
    ratio = xrange/yrange if xrange/yrange > 1 else yrange/xrange
    ratio = round(ratio)
    # Calculate the PES on a very dense and uniform grid. Useful for further 
    # analysis (MEP, shear strength) and to plot the PES
    coordinates = generate_uniform_grid(cell*ratio, density=point_density)
    E_new = rbf(coordinates[:, 0], coordinates[:, 1])
    # E_new = rbf(coordinates[:,:2])
    pes_data = np.column_stack([coordinates[:, :2], E_new])
    
    min_x = min(coordinates[:,0])
    max_x = max(coordinates[:,0])
    len_x = max_x - min_x
    dist_x = len_x/(point_density*10)
    min_y = min(coordinates[:,1])
    max_y = max(coordinates[:,1])
    len_y = max_y - min_y
    dist_y = len_y/(point_density*10)
    grid_x = np.arange(min_x, max_x, dist_x)
    grid_y = np.arange(min_y, max_y, dist_y)
    X, Y = np.meshgrid(grid_x, grid_y)
    
    # xy = np.vstack([X.ravel(), Y.ravel()]).T
    # Z = rbf(xy)
    # Z = np.reshape(Z, X.shape, order='C')
    
    Z = rbf(X, Y)
    
    to_plot = [X, Y, Z]
    
    return rbf, E_list, pes_data, data, to_plot


# =============================================================================
# UTILITY FOR THE PES
# =============================================================================

def remove_duplicates(data, rounding_decimal=5):
        
    xy = np.round(data[:, :2], decimals=rounding_decimal)
    E_dict = {}
    for i in xy:
        duplicates = np.where((xy == i).all(axis=1))[0]
        E_list = []
        for j in duplicates:
            E_list.append(data[j,2])
        E_dict[str(i)]=np.mean(E_list)
    xy_unique = np.unique(xy, axis=0)
    E=[]
    for i in xy_unique:
        E.append([E_dict[str(i)]])
    
    return np.hstack((xy_unique, np.array(E)))

def unfold_pes(hs_all, E_unique):
    """
    Unfold the energies calculated for the unique HS points of an interface,
    associating them to the replicated HS points covering the whole surface
    cell. hs_all is a dictionary, E_unique a list.

    Parameters
    ----------
        
    hs_all : dict
        Surfacial HS sites that has been unfolded (replicated) across the whole
        lattice cell of slab. 
        Ex. To the key 'ontop_1 + bridge_1' will correspond n points, spanning
        the entire plane axb of the lattice cell. Data is a (n, 3) numpy array.
    
    E : list
        Contains the energy calculated for each unique interfacial HS site.
        The energy are calculated by means of ab initio simulations (VASP).
        E must have the following structure:
            
            [ [label_1, x_1, y_1, E_1], 
              [label_2, x_2, y_2, E_2], 
              ...                      ]
        
        Ex. label should corresponds to the keys in hs_all, associated to a 
        certain shit between the lower and upper slab, e.g. 'ontop_1+bridge_1'.
        
    Returns
    -------
    E_list : list
        It's basically the same as E_unique but contains all the HS points 
        replicated on the whole interface. The structure of the list is:
            
            [ [label_1, x_1, y_1, E_1], 
              [label_2, x_2, y_2, E_2], 
              ...                      ]
        
    E_array : np.ndarray
        Numpy matrix containing the coordinates and the energy useful to 
        interpolate the PES. It's E_list without labels and with array type. 
        The structure of the matrix is:
            
            np.array([ [x_1, y_1, E_1], 
                       [x_2, y_2, E_2], 
                       ...             ])

    """

    # Initialize lists for the result
    E_list = []
    E_array = []
    
    # Extract the element
    for element in E_unique:
       label  = element[0]
       energy = element[3]
       
       # Associate each Energy to all the corresponding HS values
       for row in hs_all[label]:
          x_shift = row[0]
          y_shift = row[1]
          
          E_list.append([label, x_shift, y_shift, energy])
          E_array.append([x_shift, y_shift, energy])
          
    E_array = np.array(E_array)      
    
    return E_list, E_array

def plot_pes(data, lattice, to_fig=None):
    """
    Plot the PES and eventually save it

    """
    
    import matplotlib.pyplot as plt
    
    a = lattice[0]
    b = lattice[1]
    x = data[:, 0]
    y = data[:, 1]
    E = data[:, 2]
    
    fact=1.
    level= 43
    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    anglerot='vertical'
    shrin=1.
    zt1=plt.contourf(x, y, E, level, extent=(-fact*a, fact*a, -fact*b, fact*b), cmap=plt.cm.RdYlBu_r)
    cbar1=plt.colorbar(zt1,ax=ax,orientation=anglerot,shrink=shrin)
    cbar1.set_label(r'$E_{adh} (J/m^2)$', rotation=270, labelpad=20,fontsize=15,family='serif')
    
    ax.quiver(0. , 0., 1., 0.,scale=1.,scale_units='inches',width=0.01,color='white')
    ax.quiver(0. , 0., 0., 1.,scale=1.,scale_units='inches',width=0.01,color='white')
    ax.plot(0.,0.,'w.',ms=7)
    ax.text(0.5,-0.5,'[1 0 1]',rotation='horizontal',color='white', fontsize=14)
    ax.text(-0.5,1.,'[1 2 1]',rotation='vertical',color='white', fontsize=14)
    ax.axis([-fact*a, fact*a, -fact*b, fact*b])
    plt.xlabel(r"distance ($\AA$)",fontsize=12,family='serif')
    plt.ylabel(r"distance ($\AA$)",fontsize=12,family='serif')

    for zt1 in zt1.collections:
       zt1.set_edgecolor("face")
       zt1.set_linewidth(0.000000000001)
    
    if to_fig != None:
        plt.title("PES for " + str(to_fig), fontsize=18, family='serif')
        plt.savefig('PES_' + str(to_fig) + '.pdf', dpi=300)


# =============================================================================
# TESTING
# =============================================================================


if __name__ == '__main__':
    print('Testing the creation of a uniform grid for the PES\n')
    vectors = np.array([[3, 0, 0], [0.8, 4, 0.3]])
    a = generate_uniform_grid(vectors, density=1, pts_a=5, to_plot=True)
