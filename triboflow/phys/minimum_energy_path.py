#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:40:41 2020

Python functions to get the Minimum Energy Path (MEP) of a PES of an interface

@author: gl
"""

import numpy as np
from triboflow.phys.shear_strength import TakeDerivative, GetShearStrength_xy


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
