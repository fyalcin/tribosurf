#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 18:12:33 2020

Calculate PES, MEP, shear strength

@author: gl
"""

import numpy as np
from utility_functions import *


# =============================================================================
# MAIN
# =============================================================================


def StaticTribo(interface, hs, E):  
    """
    Main function to calculate the PES, MEP, and shear strength of an interface

    Parameters
    ----------
    interface : pymatgen.core.surface.Slab (pymatgen.core.structure.Structure)
        The interface object of which you want to calculate the PES. It is 
        needed to extract the lattice parameter used to unfold the energies.
        type(slab) could be either a Slab object or a Structure object.
        
    hs : dict
        Unfolded HS points of the interface, covering all the surface.        

    E : dict
        Contains the energy calculated for each unique surfacial HS site.

    Returns
    -------
    pes_info : list
        Relevant information about PES
        
    mep_info : list
        Relevant information about MEP
        
    ss_info : list
        Relevant information about Shear Strength
    """
    
    lattice = interface.lattice.matrix
    
    # Get the PES
    rbf, pes_dict, pes_data = GetPES(lattice, hs, E)
    data_ortho, lattice_ortho = Orthorombize(lattice, pes_data)
    
    # Get the MEP on the potential energy surface starting from a guess
    bsmep, ss_bsmep, theta = GetBSMEP(lattice_ortho, rbf)
    mep, mep_convergency = GetMEP(lattice_ortho, rbf, theta)
    
    # Calculate the MEP along the x, y, and MEP directions
    data_ss_xy, ss_xy = GetShearStrength_xy(lattice_ortho, rbf)
    data_ss_mep, ss_mep  = GetShearStrength(mep, rbf)
    
    pes_info = [pes_dict, pes_data, rbf]
    mep_info = [data_ss_mep, mep_convergency, mep, bsmep]
    ss_info = [ss_xy, ss_bsmep, ss_mep]
    
    return pes_info, mep_info, ss_info


# =============================================================================
# EVALUATION OF THE PES
# =============================================================================


def GetPES(lattice, hs_all, E, to_fig=None):
    """
    Main function to get the Potential Energy Surface (PES) for an interface. 
    The points are replicated to span a 3x3 lattice cell and are interpolated
    by using Radial Basis Functions (cubic function).
    
    ----------
    lattice : numpy.ndarray
        Vectors of the lattice cell
        
    hs : dict
        Unfolded HS points of the interface, covering all the surface.
        
    E : dict
        Contains the energy calculated for each unique surfacial HS site.
        
    to_fig : string, optional
        Name of the image that you want to save, it will be: 'name_fig'+'.pdf' 
        Suggested name: name = 'PES_' + 'Name of the interface'.
        The default is None and no image is saved.

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
    
    from scipy.interpolate import Rbf
    
    # Unfold the PES points
    data, pes_dict = UnfoldPES(E, hs_all)   
    
    # Interpolate the data with Radial Basis Function
    data_rep = ReplicatePESPoints(lattice, data, replicate_of=(3, 3) )
    rbf = Rbf(data_rep[:, 0], data_rep[:, 1], data_rep[:, 2], function='cubic')
    
    # Calculate the PES on a very dense and uniform grid. Useful for further 
    # analysis (MEP, shear strength) and to plot the PES
    coordinates = GenerateGridForPES(lattice, density=10)
    E_new = rbf(coordinates[:, :2])
    pes_data = np.column_stack([coordinates[:, :2], E_new])
    
    return rbf, pes_dict, pes_data


# =============================================================================
# UTILITY FOR THE PES
# =============================================================================


def UnfoldPES(E, hs_all):
    """
    Unfold the energies calculated for the unique HS points of an interface,
    associating them to the replicated HS points covering the whole surface
    cell.    

    Parameters
    ----------
    E : dict
        Contains the energy calculated for each Hs site of the interface.
        Ex. E may contains the key 'ontop_1 + bridge_1', corresponding to a 
        certain shift between lower and upper slab. The value of the dictionary
        is the ab initio, equilibrium energy obtained for that configuration.
        
    hs_all : dict
        Surfacial HS sites that has been unfolded (replicated) across the whole
        lattice cell of slab. 
        Ex. To the key 'ontop_1 + bridge_1' will correspond n points, spanning
        the entire plane axb of the lattice cell. Data is a (n, 3) numpy array.
        
    Returns
    -------
    pes_data : np.ndarray
        Numpy matrix containing the coordinates and the energy useful to 
        interpolate the PES. The structure of the matrix is:
            
            x[0]  y[0]  E[0]
            x[1]  y[1]  E[1]
             .     .     .
             .     .     .
             .     .     .
        
    pes_dict : dict
        Dictionary containing the coordinates and the corresponding energies
        associated to each HS point of the interface.

    """
    
    pes_dict = {}
    pes_data = []
    
    # WARNING: The elements of hs_all should not have the z coordinates.
    # Call RemoveZCoords() before using this function
    for k in hs_all.keys():
        data = np.column_stack((hs_all[k], np.full(hs_all[k], E[k])))
        pes_dict[k] = data
        pes_data.append(data)
    
    pes_data = np.concatenate(pes_data)
    
    return pes_data, pes_dict


def ReplicatePESPoints(lattice, pes_data, replicate_of=(1, 1)):
    """ 
    Replicate the PES points to cover a (n,m)-size lattice cell
    
    """
    
    n = int(replicate_of[0])
    m = int(replicate_of[1])
    
    # Check wether the number inserted are correct
    if n<=0: n=1
    if m<=0: m=1
    
    if n == 1 and m == 1:
        return pes_data
    
    else:        
        a = lattice[0, :]
        b = lattice[1, :]
        
        x = pes_data[:, 0]
        y = pes_data[:, 1]
        E = pes_data[:, 2]     
        
        x_new = np.array([])
        y_new = np.array([])
        E_new = np.array([])    
        
        for i in range(n):
                for j in range(m):
                    
                    # Replicate the x- and y- coordinates
                    x_add = x + a[0]*i + b[0]*j
                    y_add = y + a[1]*i + b[1]*j
                    
                    # Collect coordinates and energies
                    x_new = np.append(x_new, [x_new, x_add])
                    y_new = np.append(y_new, [y_new, y_add])
                    E_new = np.append(E_new, E)
        
        coordinates_new = np.column_stack([x_new, y_new, E_new])
    
        return coordinates_new


def GenerateGridForPES(lattice, density=1, pts_a=None, to_plot=False):
    """
    Generate a 2D-uniform grid of points of density=density on a lattice plane
    given by lattice[0,:]Xlattice[1,:]

    Parameters
    ----------
    lattice : numpy.ndarray
        Vectors of the lattice cell. A uniform grid of points is generated on
        the surface spanned by the first and second vector, i.e. axb.
        lattice shape is (2, 3) or (3, 3), the third vector is not necessary.
        lattice is in Angstrom units.

    density : int, optional
        Density of the grid of points that will cover the planar surface of 
        the lattice cell. Units: number of points per unit Angstrom^2
        
    pts_a : int, optional
        If this value is provided, the grid will contain pts_a points along 
        the first vector and (b/a)*pts_a along the second vector. 
        a,b : lengths of the planar lattice vectors. The default is None.
                
    to_plot : bool, optional
        Wether to display the grid of points inside the lattice cell. 
        Plot is redirected to standard output. The default is False.

    Returns
    -------
    matrix : numpy.ndarray
        Grid of points spanning the entire lattice plane.
        Format:
            x[0]  y[0]  z[0]
            y[1]  y[1]  z[1]
             .     .     .
             .     .     .
             .     .     .

    """
        
    a = lattice[0, :]
    b = lattice[1, :]
    a_mod = np.sqrt(a[0]**2. + a[1]**2. + a[2]**2.)
    b_mod = np.sqrt(b[0]**2. + b[1]**2. + b[2]**2.)
    ratio = b_mod/a_mod
    
    # Calculate the number of points for each lattice vector
    if pts_a == None:
        N_tot = round(density * a_mod * b_mod)
        n_a = int(round( np.sqrt( N_tot/ratio )))
        n_b = int(round( ratio*n_a ))
    else:
        n_a = pts_a
        n_b = int(round( ratio*n_a ))
    
    # Obtain the displacements along a and b
    dist_a_x = a[0]/n_a 
    dist_a_y = a[1]/n_a
    dist_a_z = a[2]/n_a
    dist_b_x = b[0]/n_b
    dist_b_y = b[1]/n_b
    dist_b_z = b[2]/n_b
    
    # Create the grid
    matrix = np.zeros((n_a*n_b, 3))
    k = 0
    for i in range(0, n_a):
        for j in range(0, n_b):
            matrix[k, 0] = i*dist_a_x + j*dist_b_x
            matrix[k, 1] = i*dist_a_y + j*dist_b_y
            matrix[k, 2] = i*dist_a_z + j*dist_b_z
            k += 1
    if to_plot:
        Plot_UniformGrid(lattice, matrix, n_a, n_b)

    return matrix


def Orthorombize(lattice, pes_data):
    """
    Take the replicated points of the pes and cut them in a squared shape.
    TODO : Improve the code and VECTORIZE

    """
    
    a = lattice[0, :]
    b = lattice[0, :]
    
    if np.sign(a[0]) == np.sign(b[0]):
        if a[0] > 0:
            x_up = a[0] + b [0]
            x_dw = 0
        else:
            x_up = 0
            x_dw = a[0] + b [0]           
    else:
        x_up =  max(abs(a[0]), abs(b[0]))
        x_dw =  min(abs(a[0]), abs(b[0]))
    
    if np.sign(a[1]) == np.sign(b[1]):
        if a[1] > 0:
            y_up = a[1] + b[1]
            y_dw = 0
        else:
            y_up = 0
            y_dw = a[1] 
    else:
        y_up =  max(abs(a[1]), abs(b[1]))
        y_dw =  min(abs(a[1]), abs(b[1]))
        
    index_x = pes_data[:, 0] <= 2*x_up and pes_data[:, 0] >= 2*x_dw
    index_y = pes_data[:, 1] <= 2*y_up and pes_data[:, 1] >= 2*y_dw 
    index = index_x * index_y
    
    orthorombic =  []
    for i, row in enumerate(pes_data):
        if index[i] == True:
            orthorombic.append(row)
    
    orthorombic = np.column_stack(orthorombic)
    cell = np.array([[x_up, y_dw], [x_dw, y_up]])
    
    return orthorombic, cell
    

# =============================================================================
# EVALUATION OF THE MEP
# =============================================================================


def GetMEP(lattice, rbf, theta=0., params=None):
    """
    Calculate the Minimum Energy Path on top of the Potential Surface rbf

    Parameters
    ----------
    lattice : numpy.ndarray
        Vectors of the lattice cell containing the studied interface or surface
        
    rbf : scipy.interpolate.rbf.Rbf
        Contain the information of the interpolation of the potential energy.
        
    theta : float, optional
        Slope of the starting string with respect to the x axis. The MEP will 
        be calculated starting from this directoin. To select a meaningful
        starting angle, run GetBSMEP first. Units are radiants.
        The default value is 0, i.e. x-axis
        
    params : dict, optional
        Computational parameters needed to the algorightm computing the MEP. 
        The default is None. If left to None, the following values are used:
        
        n = 101
        fac = [0.75, 0.75]
        h = 0.001
        nstepmax = 99999
        tol1 = 1e-7
        delta = 0.01
        
        n : Number of points along the string.
        fac : the string extension is (-fac, fac)*v, v=x,y. fac = [facx, facy] 
        h : time-step, limited by ODE step but independent from n1 
        nstepmax : max number of iteration to get the MEP to converge
        tol1 : tolerance for the convergence of the MEP
        delta : discretized step along x and y for integration 
        
        You can change these values by providing a dictionary containing where
        the name of the variables are the keys. Ex. {'n' : 200}
        You can change all or just some of the parameters.

    Returns
    -------
    mep : numpy.dnarray
        Set of coordinates along the MEP [x, y]. Run rbf(mep) to see the 
        potential energy profile along the MEP.
        
    mep_convergency : tuple
        Contains information about the algorithm convergence, i.e. (nstep,tol).
        nstep : number of steps done by the algorithm
        tol : final difference between the points of the string

    """

    from scipy.interpolate import interp1d, Rbf
    
    # WARNING: A squared lattice need to be provided
    a = lattice[0, :]
    b = lattice[1, :]
    alat_x = a[0]
    alat_y = b[1]
    
    # Initialize the parameters to run the algorithm 
    n = 101
    fac = [0.75, 0.75]
    h = 0.001
    nstepmax = 99999
    tol1 = 1e-7
    delta = 0.01
    
    # Check wether some parameters are inserted by the user
    if params != None and isinstance(params, dict):
        for k in params.keys():
            if k == 'n':
                n = params[k]
            elif k == 'fac':
                fac = params[k]
            elif k == 'h':
                h = params[k]
            elif k == 'nstepmax':
                nstepmax = params[k]
            elif k == 'tol':
                tol1 = params[k]
            elif k == 'delta':
                delta = params[k]    
    
    facx=fac[0]
    facy=fac[1]
    xa =-facx*alat_x
    ya =-facy*alat_y
    xb = facx*alat_x
    yb = facy*alat_y    
    g = np.linspace(0,1,n)

    if theta == np.pi/2 or theta == 3*np.pi/2: # y direction
        x = np.zeros(n)
        y = (yb-ya)*g+ya
    else: # best starting direction or x direction
        x = (xb-xa)*g+xa
        y = np.tan(theta)*x.copy()
    
    dx = x - np.roll(x,1)
    dy = y - np.roll(y,1)
    dx[0]=0.
    dy[0]=0.
    lxy  = np.cumsum(np.sqrt(dx**2+dy**2))
    lxy /= lxy[n-1]
    xf = interp1d(lxy,x,kind='cubic')
    x  =  xf(g)
    yf = interp1d(lxy,y,kind='cubic')
    y  =  yf(g)
    
    # Main loop
    for nstep in range(int(nstepmax)):
        # Calculation of the x and y-components of the force.
        # dVx and dVy are derivative of the potential
        x += delta
        tempValp=rbf(x,y)
        x -= 2.*delta
        tempValm=rbf(x,y)
        dVx = 0.5*(tempValp-tempValm)/delta
        x += delta
        y += delta
        tempValp=rbf(x,y)
        y -= 2.*delta
        tempValm=rbf(x,y)
        y += delta
        dVy = 0.5*(tempValp-tempValm)/delta

        x0 = x.copy()
        y0 = y.copy()
        # string steps:
        # 1. evolve
        xt = x- h*dVx
        yt = y - h*dVy
        # 2. derivative
        xt += delta
        tempValp=rbf(xt,yt)
        xt -= 2.*delta
        tempValm=rbf(xt,yt)
        dVxt = 0.5*(tempValp-tempValm)/delta
        xt += delta
        yt += delta
        tempValp=rbf(xt,yt)
        yt -= 2.*delta
        tempValm=rbf(xt,yt)
        yt += delta
        dVyt = 0.5*(tempValp-tempValm)/delta

        x -= 0.5*h*(dVx+dVxt)
        y -= 0.5*h*(dVy+dVyt)
        # 3. reparametrize  
        dx = x-np.roll(x,1)
        dy = y-np.roll(y,1)
        dx[0] = 0.
        dy[0] = 0.
        lxy  = np.cumsum(np.sqrt(dx**2+dy**2))
        lxy /= lxy[n-1]
        xf = interp1d(lxy,x,kind='cubic')
        x  =  xf(g)
        yf = interp1d(lxy,y,kind='cubic')
        y  =  yf(g)
        tol = (np.linalg.norm(x-x0)+np.linalg.norm(y-y0))/n
        if tol <= tol1:
           break
      
    mep = np.column_stack([x, y])
    mep_convergency = (nstep, tol)
    
    return mep, mep_convergency


def GetBSMEP(lattice, rbf, params=None):
    """
    Evaluate the best starting MEP from which start the calculation of the MEP.

    Parameters
    ----------
    lattice : numpy.ndarray
        Vectors of the lattice cell containing the studied interface or surface
        
    rbf : scipy.interpolate.rbf.Rbf
        Contain the information of the interpolation of the potential energy.
    
    params : dict, optional
        Computational parameters needed to the algorightm computing the MEP. 
        The default is None. If left to None, the following values are used:
        
        n = 1000
        fac = [1.5, 1.5]
        delta = 0.001
        delta_theta = 0.5
        
        n : Number of points along the string.
        fac : the string extension is (-fac, fac)*v, v=x,y. fac = [facx, facy] 
        delta : discretized step along x and y for integration 
        delta_theta : angle interval (degrees) to calculate the BSMEP
        
        You can change these values by providing a dictionary containing where
        the name of the variables are the keys. Ex. {'n' : 200}
        You can change all or just some of the parameters.

    Returns
    -------
    bsmep : TYPE
        DESCRIPTION.
    ss_bsmep : TYPE
        DESCRIPTION.
    theta : float
        Best orientation for starting calculating the MEP. Units are radiants

    """
    
    n = 1000
    fac = [1.5, 1.5]
    delta = 0.001
    delta_theta = 0.5
    
    if params != None and isinstance(params, dict):
        for k in params:
            if k == 'n':
                n = params[k]
            elif k == 'fac':
                fac = params[k]
            elif k == 'delta':
                delta = params[k]
            elif k == 'delta_theta':
                delta_theta = params[k]
    
    facx=fac[0]
    facy=fac[1]
    alat_x = lattice[0, 0]
    alat_y = lattice[1, 1]   
    xa =-facx*alat_x
    ya =-facy*alat_y
    xb = facx*alat_x
    yb = facy*alat_y    
    
    g = np.linspace(0,1,n)
    x = (xb-xa)*g+xa
    
    delta_theta *= np.pi/180 # Convert the angle to radiants
    rep = int(np.pi/delta_theta)    
    data_ss = []
    data_th = []   
    theta = 0
    i = 0
    
    ss_x, ss_y = GetShearStrength_xy(lattice, rbf)
    data_ss.append(float(ss_x))
    data_ss.append(float(ss_y))
    data_th.append(0.)
    data_th.append(np.pi/2.)
    
    while i < rep:
        theta += delta_theta
        if abs(theta*180/np.pi - 90) < 1e-15:
            i += 1 
            pass
        else:
            m = np.tan(theta)
            y = m*x.copy()                       

            zdev=np.zeros(len(x))        
            for i in range(len(x)):
                coordx=x[i]
                coordy=y[i]
                zdev[i]=TakeDerivative(rbf, coordx, coordy, m, delta)        
            #Shear strength in GPa
            ss = np.amax(np.abs(zdev))*10.0
            
            data_th.append(theta)
            data_ss.append(float(ss))
            i +=1
    
    index = data_ss.index(np.amin(np.abs(data_ss)))
    theta = data_th[index]
    ss_bsmep = data_ss[index]
    
    if theta == np.pi/2. :
        x = np.zeros(n)
        y = (yb-ya)*g+xa
    elif theta == 0. :
        x = (xb-xa)*g+xa
        y = np.zeros(n)
    else:
        x = (xb-xa)*g+xa
        y = np.tan(theta)*x.copy()
    
    # Tiny perturbation of the string
    x += np.random.rand(len(x))*alat_x/100
    y += np.random.rand(len(y))*alat_y/100
    bsmep = np.column_stack([x, y])
    
    return bsmep, ss_bsmep, theta


# =============================================================================
# EVALUATION OF THE SHEAR STRENGTH
# =============================================================================


def GetShearStrength(coords, rbf, delta=0.01):
    """
    Calculate the shear strength given a path and a potential energy surface.

    Parameters
    ----------
    coords : numpy.ndarray
        Coordinates [x, y] of the path along which you evaluate shear strength.
        
    rbf : scipy.interpolate.rbf.Rbf
        Contain the information of the interpolation of the potential energy.
        
    delta : TYPE, optional
        discretized step along x and y for integration. Tuning this value may
        vary slightly the final result. The default is 0.01.

    Returns
    -------
    data_ss_mep : numpy.ndarray
        Profile of potential energy and forces along the MEP.
    
    ss : float
        The shear strenth along the MEP.
    
    TODO : check what is stored in data_ss_mep and keep only important stuff

    """
    
    x = coords[:, 0]
    y = coords[:, 1]
    n = len(x)
    
    dx = x-np.roll(x,1)
    dy = y-np.roll(y,1)
    dx[0] = 0.
    dy[0] = 0.
    tx = 0.5*(np.roll(x, -1)-np.roll(x, 1))
    ty = 0.5*(np.roll(y, -1)-np.roll(y, 1))    
    # potential computed as integral of projection of gradV on string tangent
    Vz = np.zeros(n)
    #derivative of the potential
    x += delta
    tempValp=rbf(x,y)
    x -= 2.*delta
    tempValm=rbf(x,y)
    dVx = 0.5*(tempValp-tempValm)/delta
    x += delta
    y += delta
    tempValp=rbf(x,y)
    y -= 2.*delta
    tempValm=rbf(x,y)
    y += delta
    dVy = 0.5*(tempValp-tempValm)/delta

    tforce= -(tx*dVx+ty*dVy)
    force= tforce/np.sqrt(tx**2+ty**2)
    
    for i in range(n-1):
        Vz[i+1]=Vz[i] - 0.5*(tforce[i]+tforce[i+1])
        
    Vz -= np.min(Vz)
    Ve = rbf(x,y)
    Ve -= np.min(Ve)
    lxy  = np.cumsum(np.sqrt(dx**2+dy**2))
    data_ss_mep =np.stack((lxy, dVx, dVy, Vz, Ve, force), axis=-1)
    
    ss_min = 10.*np.min(force)
    ss_max = 10.*np.max(force)
    ss_mep = np.max(abs(ss_min), abs(ss_max))
    
    return data_ss_mep, ss_mep
        

def GetShearStrength_xy(lattice, rbf, params=None):   
    """
    Calculate the shear strength along the x and y directions of the cell.
    Simplified version of GetShearStrength.
    
    TODO : generalize the function in order to calculate the SS along any 
    straigth line
    
    """
    
    delta = 0.01
    npoints = 300
    
    if params != None and isinstance(params, dict):
        for k in params:
            if k == 'delta':
                delta = params[k]
            elif k == 'npoints':
                npoints = params[k]   
    
    alat_x = lattice[0, 0]
    alat_y = lattice[1, 1]
    
    x=np.arange(-1.5*alat_x, 1.5*alat_x, alat_x/npoints)
    y=np.arange(-1.5*alat_y, 1.5*alat_y, alat_y/npoints)   
    zdev_x=np.zeros(len(x))
    zdev_y=np.zeros(len(y))
    
#    zdev_x, zdev_y = derivative_xy(rbf, x, y, delta)   
    for i in range(len(x)):
        coordx=x[i]
        coordy=y[i]
        zdev_x[i] = TakeDerivative(rbf, coordx, coordy, m=0, delta=delta)
        zdev_y[i] = TakeDerivative(rbf, coordx, coordy, m=None, delta=delta)
    
    #Shear strength in GPa    
    ss_x = np.amax(np.abs(zdev_x))*10.0
    ss_y = np.amax(np.abs(zdev_y))*10.0
    ss_xy = (ss_x, ss_y)
    
    data_ss_xy = np.stack((zdev_x, zdev_y), axis=-1)
    
    return data_ss_xy, ss_xy


# =============================================================================
# UTILITY FOR THE SHEAR STRENGTH
# =============================================================================


def TakeDerivative(rbf, coordx, coordy, m=None ,delta=0.01):
    """
    Inner function used to calculate the shear strength

    """
    
    if m == None: # Derive along y
        coordx_1 = coordx
        coordx_2 = coordx
        coordy_1 = coordy - delta
        coordy_2 = coordy + delta
    elif m == 0: # Derive along x
        coordx_1 = coordx - delta
        coordx_2 = coordx + delta
        coordy_1 = coordy
        coordy_2 = coordy
    else: # Derive along the straight line with slope m
        coordx_1 = coordx - m * delta 
        coordx_2 = coordx + m * delta
        coordy_1 = coordy
        coordy_2 = coordy
    
    # Calculate the derivative
    V_1 = rbf(coordx_1, coordy_1)
    V_2 = rbf(coordx_2, coordy_2)
    zdev = 0.5 * (V_2 - V_1) / delta
    
    return zdev


# OLD VERSIONS OF THE DERIVATIVE. IT WORKED

# def derivative_xy(rbf, coordx, coordy, delta):
#     coordx=coordx+delta
#     tempValp=rbf(coordx,0.)
#     coordx=coordx - 2.*delta
#     tempValm=rbf(coordx,0.)
#     zdev_x=0.5*(tempValp-tempValm)/delta
    
#     coordy=coordy+delta
#     tempValp=rbf(0.,coordy)
#     coordy=coordy - 2.*delta
#     tempValm=rbf(0.,coordy)
#     zdev_y=0.5*(tempValp-tempValm)/delta
    
#     return zdev_x, zdev_y


# def derivative_bsmep(rbf, coordx, coordy, m, delta):
#     coordx=coordx+delta
#     coordy=coordy+m*delta
#     tempValp=rbf(coordx,coordy)
    
#     coordx=coordx - 2.*delta
#     coordy=coordy - 2.*delta*m
#     tempValm=rbf(coordx,coordy)
#     zdev=0.5*(tempValp-tempValm)/delta
    
#     return zdev


# =============================================================================
# TESTING
# =============================================================================


if __name__ == '__main__':
    print('Testing the creation of a uniform grid for the PES\n')
    vectors = np.array([[3, 0, 0], [0.8, 4, 0.3]])
    a= GenerateGridForPES(vectors, density=1, pts_a=5, to_plot=True)