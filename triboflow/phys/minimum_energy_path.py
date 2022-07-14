#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:40:41 2020

Python functions to get the Minimum Energy Path (MEP) of a PES of an interface

The module contains the following functions:

    - get_mep
    - get_bs_mep
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
from scipy.interpolate import interp1d

from triboflow.phys.shear_strength import take_derivative, get_shear_strength_xy


# =============================================================================
# EVALUATION OF THE MEP
# =============================================================================

def get_initial_strings(extended_energy_list, xlim, ylim, npts=100, add_noise=0.01):
    energies = extended_energy_list.copy()
    energies[:, 2] = energies[:, 2] - min(energies[:, 2])

    minima = [x[:2] for x in energies if (x[2] == 0.0 and
                                          x[0] >= 0.1 and
                                          x[0] < xlim - 0.1 and
                                          x[1] >= 0.1 and
                                          x[1] < ylim - 0.1)]

    start = [xlim, ylim]
    end = [0.0, 0.0]
    end_x = [0.0, ylim]
    end_y = [0.0, 0.0]

    # make a path mostly diagonal through the plotting region
    for p in minima:
        if np.linalg.norm(p) < np.linalg.norm(start):
            start = p
        if np.linalg.norm(p) > np.linalg.norm(end):
            end = p

    # make a path mostly in x direction through the plotting region
    for p in minima:
        if p[0] > end_x[0] and np.isclose(p[1], start[1], atol=0.01):
            end_x = p
    # make a path mostly in y direction through the plotting region
    for p in minima:
        if p[1] > end_y[1] and np.isclose(p[0], start[0], atol=0.01):
            end_y = p

    string_d = np.linspace(start, end, npts)
    string_x = np.linspace(start, end_x, npts)
    string_y = np.linspace(start, end_y, npts)
    if add_noise:
        noise = np.random.normal(0, add_noise, string_d.shape)
        string_d += noise
        string_x += noise
        string_y += noise
    return string_d, string_x, string_y


def new_evolve_string(string, rbf, nstepmax=99999, mintol=1e-8, delta=0.005, h=0.01):
    """
    Find a minumum energy path from an initial string.
    
    Simplified zero temperature string method as described by Weinan E et al.
    J. Chem. Phys. 126, 164103 (2007)
    

    Parameters
    ----------
    string : numpy.ndarray
        The initial string
    rbf : scipy.interpolate._rbfinterp.RBFInterpolator
        Radial Basis Function that describes the potential in which the string
        is evolved
    nstepmax : int, optional
        Maximum number of time steps. The default is 99999.
    mintol : float, optional
        Minimal tolerance to define convergence. The default is 1e-8.
    delta : float, optional
        Delta parameter for the gradient calculation. The default is 0.005.
    h : float, optional
        timestep. The default is 0.01.

    Returns
    -------
    dict
        Minimum energy path as a numpy array as well as the convergence
        parameters.

    """
    n = len(string)
    g = np.linspace(0, 1, n)
    x = string[:,0]
    y = string[:,1]
    
    for nstep in range(int(nstepmax)):
        xp = x + delta
        xm = x - delta
        yp = y + delta
        ym = y - delta
        
        #energy = rbf(string)
        energy_xm = rbf(np.stack((xm, y), axis=1))
        energy_xp = rbf(np.stack((xp, y), axis=1))
        energy_ym = rbf(np.stack((x, ym), axis=1))
        energy_yp = rbf(np.stack((x, yp), axis=1))
        
        # pot_x = np.stack((energy_xm, energy, energy_xp), axis=1)
        # pot_y = np.stack((energy_ym, energy, energy_yp), axis=1)
        pot_x = np.stack((energy_xm, energy_xp), axis=1)
        pot_y = np.stack((energy_ym, energy_yp), axis=1)
        #cannot get this to work in one step
        #gradient = np.gradient(energy, string[:,0], string[:,1], axis=0)
        # gradx = np.gradient(pot_x, delta, axis=1)[:,1]
        # grady = np.gradient(pot_y, delta, axis=1)[:,1]
        gradx = np.gradient(pot_x, 2*delta, axis=1)[:,0]
        grady = np.gradient(pot_y, 2*delta, axis=1)[:,0]
        
        string_old = string.copy()
        #evolve the string
        string_new = np.stack((string[:,0] - h * gradx, string[:,1] - h * grady), axis=1)
        
        
        # 3. reparametrize  
        x = string_new[:,0]
        y = string_new[:,1]
        dx = np.ediff1d(x, to_begin=0)
        dy = np.ediff1d(y, to_begin=0)
        lxy = np.cumsum(np.sqrt(dx ** 2 + dy ** 2))
        lxy /= lxy[n - 1]
        xf = interp1d(lxy, x, kind='cubic')
        x = xf(g)
        yf = interp1d(lxy, y, kind='cubic')
        y = yf(g)
        string = np.stack((x, y), axis=1)
        
        #check for convergence
        tol = (np.linalg.norm(string - string_old)) / n
        if tol <= mintol:
            break
        
    mep_convergency = (nstep, tol)
    return {'mep': string, 'convergence': mep_convergency}
    
def evolve_string(string, rbf, nstepmax=99999, mintol=1e-7, delta=0.005, h=0.001):
    x = string[:, 0]
    y = string[:, 1]
    n = len(x)
    g = np.linspace(0, 1, n)
    dx = x - np.roll(x, 1)
    dy = y - np.roll(y, 1)
    dx[0] = 0.
    dy[0] = 0.
    lxy = np.cumsum(np.sqrt(dx ** 2 + dy ** 2))
    lxy /= lxy[n - 1]
    # xf = interp1d(lxy,x,kind='cubic')
    # x  =  xf(g)
    # yf = interp1d(lxy,y,kind='cubic')
    # y  =  yf(g)
    for nstep in range(int(nstepmax)):
        # Calculation of the x and y-components of the force.
        # dVx and dVy are derivative of the potential
        x += delta
        tempValp = rbf(np.stack((x, y), axis=1))
        x -= 2. * delta
        tempValm = rbf(np.stack((x, y), axis=1))
        dVx = 0.5 * (tempValp - tempValm) / delta
        x += delta
        y += delta
        tempValp = rbf(np.stack((x, y), axis=1))
        y -= 2. * delta
        tempValm = rbf(np.stack((x, y), axis=1))
        y += delta
        dVy = 0.5 * (tempValp - tempValm) / delta

        x0 = x.copy()
        y0 = y.copy()
        # string steps:
        # 1. evolve
        xt = x - h * dVx
        yt = y - h * dVy
        # 2. derivative
        xt += delta
        tempValp = rbf(np.stack((xt, yt), axis=1))
        xt -= 2. * delta
        tempValm = rbf(np.stack((xt, yt), axis=1))
        dVxt = 0.5 * (tempValp - tempValm) / delta
        xt += delta
        yt += delta
        tempValp = rbf(np.stack((xt, yt), axis=1))
        yt -= 2. * delta
        tempValm = rbf(np.stack((xt, yt), axis=1))
        yt += delta
        dVyt = 0.5 * (tempValp - tempValm) / delta

        x -= 0.5 * h * (dVx + dVxt)
        y -= 0.5 * h * (dVy + dVyt)
        # 3. reparametrize  
        dx = x - np.roll(x, 1)
        dy = y - np.roll(y, 1)
        dx[0] = 0.
        dy[0] = 0.
        lxy = np.cumsum(np.sqrt(dx ** 2 + dy ** 2))
        lxy /= lxy[n - 1]
        xf = interp1d(lxy, x, kind='cubic')
        x = xf(g)
        yf = interp1d(lxy, y, kind='cubic')
        y = yf(g)
        tol = (np.linalg.norm(x - x0) + np.linalg.norm(y - y0)) / n
        if tol <= mintol:
            break

    mep = np.column_stack([x, y])
    mep_convergency = (nstep, tol)
    return {'mep': mep, 'convergence': mep_convergency}


def get_mep(lattice, rbf, theta=0., params=None):
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
        starting angle, run get_bs_mep first. Units are radiants.
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

    from scipy.interpolate import interp1d

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

    facx = fac[0]
    facy = fac[1]
    xa = -facx * alat_x
    ya = -facy * alat_y
    xb = facx * alat_x
    yb = facy * alat_y
    g = np.linspace(0, 1, n)

    if theta == np.pi / 2 or theta == 3 * np.pi / 2:  # y direction
        x = np.zeros(n)
        y = (yb - ya) * g + ya
    else:  # best starting direction or x direction
        x = (xb - xa) * g + xa
        y = np.tan(theta) * x.copy()

    dx = x - np.roll(x, 1)
    dy = y - np.roll(y, 1)
    dx[0] = 0.
    dy[0] = 0.
    lxy = np.cumsum(np.sqrt(dx ** 2 + dy ** 2))
    lxy /= lxy[n - 1]
    xf = interp1d(lxy, x, kind='cubic')
    x = xf(g)
    yf = interp1d(lxy, y, kind='cubic')
    y = yf(g)

    # Main loop
    for nstep in range(int(nstepmax)):
        # Calculation of the x and y-components of the force.
        # dVx and dVy are derivative of the potential
        x += delta
        tempValp = rbf(x, y)
        x -= 2. * delta
        tempValm = rbf(x, y)
        dVx = 0.5 * (tempValp - tempValm) / delta
        x += delta
        y += delta
        tempValp = rbf(x, y)
        y -= 2. * delta
        tempValm = rbf(x, y)
        y += delta
        dVy = 0.5 * (tempValp - tempValm) / delta

        x0 = x.copy()
        y0 = y.copy()
        # string steps:
        # 1. evolve
        xt = x - h * dVx
        yt = y - h * dVy
        # 2. derivative
        xt += delta
        tempValp = rbf(xt, yt)
        xt -= 2. * delta
        tempValm = rbf(xt, yt)
        dVxt = 0.5 * (tempValp - tempValm) / delta
        xt += delta
        yt += delta
        tempValp = rbf(xt, yt)
        yt -= 2. * delta
        tempValm = rbf(xt, yt)
        yt += delta
        dVyt = 0.5 * (tempValp - tempValm) / delta

        x -= 0.5 * h * (dVx + dVxt)
        y -= 0.5 * h * (dVy + dVyt)
        # 3. reparametrize  
        dx = x - np.roll(x, 1)
        dy = y - np.roll(y, 1)
        dx[0] = 0.
        dy[0] = 0.
        lxy = np.cumsum(np.sqrt(dx ** 2 + dy ** 2))
        lxy /= lxy[n - 1]
        xf = interp1d(lxy, x, kind='cubic')
        x = xf(g)
        yf = interp1d(lxy, y, kind='cubic')
        y = yf(g)
        tol = (np.linalg.norm(x - x0) + np.linalg.norm(y - y0)) / n
        if tol <= tol1:
            break

    mep = np.column_stack([x, y])
    mep_convergency = (nstep, tol)

    return mep, mep_convergency


def get_bs_mep(lattice, rbf, params=None):
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

    facx = fac[0]
    facy = fac[1]
    alat_x = lattice[0, 0]
    alat_y = lattice[1, 1]
    xa = -facx * alat_x
    ya = -facy * alat_y
    xb = facx * alat_x
    yb = facy * alat_y

    g = np.linspace(0, 1, n)
    x = (xb - xa) * g + xa

    delta_theta *= np.pi / 180  # Convert the angle to radiants
    rep = int(np.pi / delta_theta)
    data_ss = []
    data_th = []
    theta = 0
    i = 0

    ss_x, ss_y = get_shear_strength_xy(lattice, rbf)
    data_ss.append(float(ss_x))
    data_ss.append(float(ss_y))
    data_th.append(0.)
    data_th.append(np.pi / 2.)

    while i < rep:
        theta += delta_theta
        if abs(theta * 180 / np.pi - 90) < 1e-15:
            i += 1
            pass
        else:
            m = np.tan(theta)
            y = m * x.copy()

            zdev = np.zeros(len(x))
            for i in range(len(x)):
                coordx = x[i]
                coordy = y[i]
                zdev[i] = take_derivative(rbf, coordx, coordy, m, delta)
                # Shear strength in GPa
            ss = np.amax(np.abs(zdev)) * 10.0

            data_th.append(theta)
            data_ss.append(float(ss))
            i += 1

    index = data_ss.index(np.amin(np.abs(data_ss)))
    theta = data_th[index]
    ss_bsmep = data_ss[index]

    if theta == np.pi / 2.:
        x = np.zeros(n)
        y = (yb - ya) * g + xa
    elif theta == 0.:
        x = (xb - xa) * g + xa
        y = np.zeros(n)
    else:
        x = (xb - xa) * g + xa
        y = np.tan(theta) * x.copy()

    # Tiny perturbation of the string
    x += np.random.rand(len(x)) * alat_x / 100
    y += np.random.rand(len(y)) * alat_y / 100
    bsmep = np.column_stack([x, y])

    return bsmep, ss_bsmep, theta
