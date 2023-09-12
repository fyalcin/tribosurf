#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Michael Wolloch"
__copyright__ = "Copyright 2022, M. Wolloch, HIT, FWF, University of Vienna"
__credits__ = "Code readapted from our past work by Gabriele Losi and Mauro Ferrario, MIT license, https://github.com/mcrighi/interface-workflow,"
__contact__ = "michael.wolloch@univie.ac.at"
__date__ = "September 16th, 2022"

import math as m
import multiprocessing
import numpy as np
from mep.models import Model
from mep.neb import NEB
from mep.path import Image, Path
from scipy.interpolate import interp1d


# =============================================================================
# EVALUATION OF THE MEP
# =============================================================================


def evolve_mep(
    string_dict, rbf, method="neb", max_iter=99999, neb_forcetol=1e-3
):
    """
    Compute minimum energy paths for all strings in the sting_dictionary.

    Either use the nudged elastiv band ('neb', default) or zero temperature
    string method ('zts') to find the MEPs. All three directions are
    optimized in parallel, but depending on the PES, this might still take
    around a couple of minutes or half an hour...

    Returns
    -------
    None.

    """

    if method == "zts":
        # evolve the stings parallel by using a Pool
        pool_inputs = []
        mep_list = []
        nsteps = max_iter
        for name, s in string_dict.items():
            if len(s) > 10:
                pool_inputs.append([s, rbf, nsteps])
                mep_list.append(name)
        if len(pool_inputs) > 0:
            try:
                evolve_pool = multiprocessing.Pool(len(pool_inputs))
                meps = evolve_pool.starmap(new_evolve_string, pool_inputs)
            except:
                print(
                    "Parallel ZTS method failed and crashed. Swiching to sequential."
                )
                # try running sequentially because pool might be running into
                # a recursion limit on some systems
                meps = []
                for inputs in pool_inputs:
                    mep = new_evolve_string(
                        string=pool_inputs[0],
                        rbf=pool_inputs[1],
                        nstepmax=pool_inputs[2],
                    )
                    meps.append(mep)
        else:
            meps = None

    elif method == "neb":
        pool_inputs = []
        mep_list = []
        nsteps = max_iter

        model = RBFModel(rbf)
        for name, s in string_dict.items():
            if len(s) > 10:
                pool_inputs.append(
                    [model, Path(s, 1), nsteps, neb_forcetol]
                )
                mep_list.append(name)
        if len(pool_inputs) > 0:
            try:
                evolve_pool = multiprocessing.Pool(len(pool_inputs))
                meps = evolve_pool.starmap(run_neb, pool_inputs)
            except:
                print(
                    "Parallel NEB method failed and crashed. Swiching to sequential."
                )
                # try running sequentially because pool might be running into
                # a recursion limit on some systems
                meps = []
                for inputs in pool_inputs:
                    mep = run_neb(
                        model=inputs[0],
                        path=inputs[1],
                        nsteps=inputs[2],
                        tol=inputs[3],
                    )
                    meps.append(mep)
        else:
            meps = None
    else:
        print(
            f'WARNING: method "{method}" is not supported. Choose either "neb" or "zts".'
        )
        return None

    mep = {}
    for i, name in enumerate(mep_list):
        mep[name] = meps[i]

    return mep


def find_minima(extended_energy_list, xlim, ylim, border_padding):
    """
    Return global minimas in an energy list.

    Parameters
    ----------
    extended_energy_list : array
        High symmetry points with x and y coordinates and the total energy
        as third column.
    xlim : float
        maximum x value for the minima search
    ylim : float
        maximum y value for the minima search
    border_padding : float, optional
        Include only minima at least border_padding away from the borders of
        the plotting region (in Angstroms). The default is 0.1

    Returns
    -------
    minima : np.array
        x and y positions of the minima of the PES.

    """
    energies = extended_energy_list.copy()
    energies[:, 2] = energies[:, 2] - min(energies[:, 2])

    minima = [
        x[:2]
        for x in energies
        if (
            x[2] == 0.0
            and x[0] >= border_padding
            and x[0] < xlim - border_padding
            and x[1] >= border_padding
            and x[1] < ylim - border_padding
        )
    ]
    return minima


def get_initial_strings(
    extended_energy_list,
    xlim,
    ylim,
    point_density=20,
    add_noise=0.01,
    border_padding=0.1,
):
    """
    Make 3 straight strings that connect minima of the PES.

    One string each is set up be parallel to the cartesian x and y directions,
    while the third one is roughly diagonal. The number of points is controlled
    by the density parameter, and noise can be added to the strings to ensure
    that some forces are present even if the string is located alongside a
    path were no normal forces are found for symmetry reasons.

    Parameters
    ----------
    extended_energy_list : numpy array
        High symmetry points with x and y coordinates and the total energy
        as third column.
    xlim : float
        maximum x value for the minima search
    ylim : float
        maximum y value for the minima search
    point_density : int, optional
        Number of points/images per Angstrom. The default is 20.
    add_noise : float or bool, optional
        Set either to False or a float value in Angstrom. The default is 0.01.
    border_padding : float, optional
        paths start only so far away from the borders of the plotting region
        (in Angstroms). The default is 0.1

    Returns
    -------
    string_d : array
        Diagonal string
    string_x : array
        String in x direction
    string_y : array
        String in y direction

    """

    minima = find_minima(extended_energy_list, xlim, ylim, border_padding)

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

    npts_x = m.ceil(
        (np.linalg.norm(end_x) - np.linalg.norm(start)) * point_density
    )
    npts_y = m.ceil(
        (np.linalg.norm(end_y) - np.linalg.norm(start)) * point_density
    )
    npts_d = m.ceil(
        (np.linalg.norm(end) - np.linalg.norm(start)) * point_density
    )

    string_d = np.linspace(start, end, npts_d)
    string_x = np.linspace(start, end_x, npts_x)
    string_y = np.linspace(start, end_y, npts_y)
    if add_noise:
        string_d += np.random.normal(0, add_noise, string_d.shape)
        string_x += np.random.normal(0, add_noise, string_x.shape)
        string_y += np.random.normal(0, add_noise, string_y.shape)
    return string_d, string_x, string_y


def numgrad(string, rbf, delta=0.002):
    """
    Numerically compute the gradient of a potential given by rbf at the points in string.

    :param string: The string in the potential
    :type string: numpy array

    :param rbf: Method to interpolate the potential
    :type rbf: scipy.interpolate._rbfinterp.RBFInterpolator or other callable

    :param delta: step size for the gradient evaluation. The default is 0.002.
    :type delta: float, optional

    :return: The gradient at each string point.
    :rtype: numpy array

    """
    x = string[:, 0]
    y = string[:, 1]

    xp = x + delta
    xm = x - delta
    yp = y + delta
    ym = y - delta

    energy = rbf(string)
    energy_xm = rbf(np.stack((xm, y), axis=1))
    energy_xp = rbf(np.stack((xp, y), axis=1))
    energy_ym = rbf(np.stack((x, ym), axis=1))
    energy_yp = rbf(np.stack((x, yp), axis=1))

    potx = np.stack((energy_xm, energy, energy_xp), axis=1)
    poty = np.stack((energy_ym, energy, energy_yp), axis=1)

    gradientx = np.gradient(potx, delta, axis=1)[:, 1]
    gradienty = np.gradient(poty, delta, axis=1)[:, 1]

    return np.stack((gradientx, gradienty), axis=1)


def reparametrize_string_with_equal_spacing(string, nr_of_points):
    """
    Reparametrize the string so that the arc length are constant again.

    Parameters
    ----------
    string : np.array
        String with non-equal arc lengths.
    nr_of_points : int
        New number of points for the string

    Returns
    -------
    np.array
        String with equal arc lengths.

    """
    g = np.linspace(0, 1, nr_of_points)
    x = string[:, 0]
    y = string[:, 1]
    dx = np.ediff1d(x, to_begin=0)
    dy = np.ediff1d(y, to_begin=0)
    lxy = np.cumsum(np.sqrt(dx**2 + dy**2))  # lxy[n-1] = sum(lxy[:n])
    lxy /= lxy[-1]  # rescale distance between points to [0,1] interval
    xf = interp1d(lxy, x, kind="cubic")  # interpolate x=f(lxy)
    x = xf(
        g
    )  # since g is evenly spaced, now the new points are evenly distributed
    yf = interp1d(lxy, y, kind="cubic")  # interpolate y=f(lxy)
    y = yf(
        g
    )  # since g is evenly spaced, now the new points are evenly distributed
    return np.stack((x, y), axis=1)


def new_evolve_string(
    string, rbf, nstepmax=9999, mintol=1e-7, delta=0.005, h=0.005
):
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
        Delta parameter for the gradient calculation. The default is 0.002.
    h : float, optional
        timestep. The default is 0.01.

    Returns
    -------
    dict
        Minimum energy path as a numpy array as well as info on convergence.

    """
    n = len(string)
    is_converged = False

    for nstep in range(int(nstepmax)):
        gradient = numgrad(string, rbf, delta)

        string_old = string.copy()
        # evolve the string
        string_new = string - h * gradient

        # 3. reparametrize the string so the arc length is equal everywhere.
        string = reparametrize_string_with_equal_spacing(string_new, n)

        # check for convergence
        tol = (np.linalg.norm(string - string_old)) / n
        if tol <= mintol:
            is_converged = True
            break

    return {"mep": string, "convergence": is_converged}


def run_neb(model, path, nsteps, tol):
    """
    Run the NEB model
    """
    neb = NEB(model, path)
    history = neb.run(n_steps=nsteps, force_tol=tol, verbose=False)
    mep = np.array([n.data.tolist()[0] for n in neb.path])
    return {"mep": mep, "convergence": neb.stop}


class RBFModel(Model):
    def __init__(self, rbf):
        self.rbf = rbf

    def predict_energy(self, image):
        if isinstance(image, Image):
            image = image.data
        image = np.atleast_2d(image)
        return self.rbf(image)
