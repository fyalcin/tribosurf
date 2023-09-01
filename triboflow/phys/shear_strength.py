#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:44:44 2020

Python functions to get the Shear Strength (SS) of an interface

The module contains the following functions:

    - get_shear_strength
    - get_shear_strength_xy
    - take_derivative

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna
    Credits: Code readapted from our past homogeneous workflow, MIT license,
    https://github.com/mcrighi/interface-workflow,
"""

__author__ = "Gabriele Losi"
__credits__ = "Code readapted from our past homogeneous workflow, MIT license, https://github.com/mcrighi/interface-workflow,"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 8th, 2021"

import numpy as np


# =============================================================================
# EVALUATION OF THE SHEAR STRENGTH
# =============================================================================


def get_shear_strength(coords, rbf, delta=0.01):
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

    dx = x - np.roll(x, 1)
    dy = y - np.roll(y, 1)
    dx[0] = 0.0
    dy[0] = 0.0
    tx = 0.5 * (np.roll(x, -1) - np.roll(x, 1))
    ty = 0.5 * (np.roll(y, -1) - np.roll(y, 1))
    # potential computed as integral of projection of gradV on string tangent
    Vz = np.zeros(n)
    # derivative of the potential
    x += delta
    tempValp = rbf(x, y)
    x -= 2.0 * delta
    tempValm = rbf(x, y)
    dVx = 0.5 * (tempValp - tempValm) / delta
    x += delta
    y += delta
    tempValp = rbf(x, y)
    y -= 2.0 * delta
    tempValm = rbf(x, y)
    y += delta
    dVy = 0.5 * (tempValp - tempValm) / delta

    tforce = -(tx * dVx + ty * dVy)
    force = tforce / np.sqrt(tx**2 + ty**2)

    for i in range(n - 1):
        Vz[i + 1] = Vz[i] - 0.5 * (tforce[i] + tforce[i + 1])

    Vz -= np.min(Vz)
    Ve = rbf(x, y)
    Ve -= np.min(Ve)
    lxy = np.cumsum(np.sqrt(dx**2 + dy**2))
    data_ss_mep = np.stack((lxy, dVx, dVy, Vz, Ve, force), axis=-1)

    ss_min = 10.0 * np.min(force)
    ss_max = 10.0 * np.max(force)
    ss_mep = np.max(abs(ss_min), abs(ss_max))

    return data_ss_mep, ss_mep


def get_shear_strength_xy(lattice, rbf, params=None):
    """
    Calculate the shear strength along the x and y directions of the cell.
    Simplified version of get_shear_strength.

    TODO : generalize the function in order to calculate the SS along any
    straigth line

    """

    delta = 0.01
    npoints = 300

    if params is not None and isinstance(params, dict):
        for k in params:
            if k == "delta":
                delta = params[k]
            elif k == "npoints":
                npoints = params[k]

    alat_x = lattice[0, 0]
    alat_y = lattice[1, 1]

    x = np.arange(-1.5 * alat_x, 1.5 * alat_x, alat_x / npoints)
    y = np.arange(-1.5 * alat_y, 1.5 * alat_y, alat_y / npoints)
    zdev_x = np.zeros(len(x))
    zdev_y = np.zeros(len(y))

    #    zdev_x, zdev_y = derivative_xy(rbf, x, y, delta)
    for i in range(len(x)):
        coordx = x[i]
        coordy = y[i]
        zdev_x[i] = take_derivative(rbf, coordx, coordy, m=0, delta=delta)
        zdev_y[i] = take_derivative(rbf, coordx, coordy, m=None, delta=delta)

    # Shear strength in GPa
    ss_x = np.amax(np.abs(zdev_x)) * 10.0
    ss_y = np.amax(np.abs(zdev_y)) * 10.0
    ss_xy = (ss_x, ss_y)

    data_ss_xy = np.stack((zdev_x, zdev_y), axis=-1)

    return data_ss_xy, ss_xy


# =============================================================================
# UTILITY FOR THE SHEAR STRENGTH
# =============================================================================


def take_derivative(rbf, coordx, coordy, m=None, delta=0.01):
    """
    Inner function used to calculate the shear strength

    """

    if m is None:  # Derive along y
        coordx_1 = coordx
        coordx_2 = coordx
        coordy_1 = coordy - delta
        coordy_2 = coordy + delta
    elif m == 0:  # Derive along x
        coordx_1 = coordx - delta
        coordx_2 = coordx + delta
        coordy_1 = coordy
        coordy_2 = coordy
    else:  # Derive along the straight line with slope m
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
