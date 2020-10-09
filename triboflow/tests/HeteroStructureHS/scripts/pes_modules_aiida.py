import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from scipy.interpolate import interp1d, Rbf, interp2d, griddata
from plot_modules import *


################################################################################
###########################      READING INPUT       ###########################
################################################################################

def read_input(inp='input_WF-shearstrength.in', out='output.txt', out2='output.txt', data='array_for_2dpes.dat', ):
    from matplotlib import use
    use('Agg')
    import numpy as np
    import os
    from aiida.orm import load_node
    
    if os.path.isfile(data):
        array_accumulo = np.genfromtxt('array_for_2dpes.dat')
    else:
        print( 'Please provide PES data as an "array_for_2dpes.dat" file\n')
        raise SystemExit
    
    try:
        inputfile = open(inp, 'r')
    except IOError:
        with open(out2, 'a') as out:
            out.write(" ")
            out.write("Error opening input file 'input_WF-shearstrength.in'. Please make sure that it is in the current directory.\n")
            out.write(" ")
        raise SystemExit
    input_dict={}
    for line in inputfile:
        line.strip()
        if not line.startswith('#') and not line.startswith('\n'):
            line = line.partition('#')[0]
            line = line.rstrip()
            name, var = line.partition("=")[::2]
            input_dict[name.strip()] = var
    inputfile.close()
    
    try:
        outfile = open(out,'r')
    except IOError:
        print(" ")
        print("Error opening output file 'output.txt'. Please make sure that it is in the current directory.\n")
        print(" ")
        raise SystemExit
    for line in outfile:
        if line.startswith('Package number of the bulk structure'):
            structure_pk = int(line.split()[-1])
            structure = load_node(structure_pk)
        if line.startswith('Selected external pressure of'):
            scale = line.split()[-1]
            scale = float(scale[:-2]) / 100.
    outfile.close()
    
    element = input_dict['element'].strip()
    lattice = input_dict['lattice'].strip()
    miller = input_dict['selected_surface'].strip()
    
    with open(out2, 'a') as out:
            out.write("Starting the PES and SS_MEP calculation for  {0}, {1}, {2}\n".format(str(element), str(lattice), str(miller)))
            
    return array_accumulo, structure, element, lattice, miller, scale

################################################################################
####################      MODULES FOR PES CALCULATION       ####################
################################################################################

def pes(array_accumulo, structure, scale=0., element='interface', lattice='', miller=''):   
    """
    Calculation of the PES
    arra_accumulo            : high simmetry points with corresponding energy
    structure                : aiida object with lattice cell information
    element, lattice, miller : information about the cell
    scale                    : rescaling factor of lattice due to pressure
    """
    
    expansion = 1.0 + scale
    x = np.array(array_accumulo[:,1])*expansion
    y = np.array(array_accumulo[:,2])*expansion
    E = np.array(array_accumulo[:,3])
    
    alat, xnew, ynew = alats_and_accumulo(array_accumulo, structure, expansion, lattice, miller)   
    rbf = Rbf(x, y, E, function='cubic')
    Enew=rbf(xnew, ynew)
    
    return alat, x, y, E, xnew, ynew, Enew, rbf

#########################

def alats_and_accumulo(array_accumulo, structure, expansion, lattice, miller, fact=1.):
    
    #Specify you own alat if needed e.g. orthorombic
    #alat_x = 3.331843818 ; alat_y = 4.444098655
    alat_x = structure.cell[0][0]
    alat_y = structure.cell[1][1]
    alat_x *= expansion
    alat_y *= expansion

    if lattice == 'bcc' and miller == '110':
        alat_x *= 2.0
        alat_y *= 2.0
    if miller == '111' or miller == '0001':
        alat_y *= 2.0
    
    #Readapt the points to the new vectors
    if miller == '100':
        xnew, ynew = np.mgrid[-fact*alat_x :fact*alat_x:260j,-fact*alat_y:fact*alat_y:260j]
    elif miller == '110':
        xnew, ynew = np.mgrid[-fact*alat_x :fact*alat_x:220j,-fact*alat_y:fact*alat_y:312j]
    elif miller == '111' or miller == '0001':
        xnew, ynew = np.mgrid[-fact*alat_x :fact*alat_x:200j,-fact*alat_y:fact*alat_y:346j]
    else:
        xnew, ynew = np.mgrid[-fact*alat_x :fact*alat_x:200j,-fact*alat_y:fact*alat_y:346j]
    
    return [alat_x, alat_y], xnew, ynew

################################################################################
################      MODULES FOR MEP AND SS CALCULATIONS       ################
################################################################################

def mep(alat, xnew, ynew, rbf, theta, n=101, fac=[0.75, 0.75], h=0.001, nstepmax=99999, tol1=1e-7, delta=0.0001):
    """
    Calculation of the MEP with the string method along a provided path
    alat        : array containing alats values [alat_x, alat_y]
    xnew, ynew  : extended x and y points on the potential energy map surface
    rbf         : potential energy surface interpolated with rbf
    theta       : starting splope of the string (with respect to x axis)
    fac         : array used for setting the string length [facx, facy]
    h           : time-step (limited by the ODE step but independent of n1)
    nstepmax    : max possible number of iterations
    tol1        : parameter used as stopping criterion
    delta       : discretized step along x and y for integration
    """
    
    facx=fac[0]
    facy=fac[1]
    alat_x = alat[0]
    alat_y = alat[1]
    
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
        # calculation of the x and y-components of the force, dVx and dVy respectively
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
       
    return x, y, nstep, tol

#########################

def shear_strength(x, y, rbf, delta=0.001):
    """
    Calculate the SS given a path and a potential energy surface
    x,y         : points of the string 
    rbf         : interpolated potential energy surface
    delta       : discretized step along x and y for integration
    """  
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
    
    SS_min = 10.*np.min(force)
    SS_max = 10.*np.max(force)
    
    return SS_max, SS_min, data_ss_mep         
        
################################################################################
###############      SS AND PEP CALCULATION ALONG X AND Y       ################
################################################################################

def shear_strength_xy(alat, rbf, npoints=300, bsmep=False):   
    
    delta=0.0001
    alat_x = alat[0]
    alat_y = alat[1]
    
    x=np.arange(-1.5*alat_x, 1.5*alat_x, alat_x/npoints)
    y=np.arange(-1.5*alat_y, 1.5*alat_y, alat_y/npoints)   
    zdev_x=np.zeros(len(x))
    zdev_y=np.zeros(len(y))
    
#    zdev_x, zdev_y = derivative_xy(rbf, x, y, delta)   
    for i in range(len(x)):
        coordx=x[i]
        coordy=y[i]
        zdev_x[i], zdev_y[i] = derivative_xy(rbf, coordx, coordy, delta)
    
    #Shear strength in GPa    
    ss_x = np.amin(zdev_x)*10.0
    ss_y = np.amin(zdev_y)*10.0
    #data=np.stack((zdev_x, zdev_y), axis=-1)

    if bsmep == True:
        ss_x = np.amax(np.abs(zdev_x))*10.0
        ss_y = np.amax(np.abs(zdev_y))*10.0
    
    return abs(ss_x), abs(ss_y)#, data

#########################

def derivative_xy(rbf, coordx, coordy, delta):
    coordx=coordx+delta
    tempValp=rbf(coordx,0.)
    coordx=coordx - 2.*delta
    tempValm=rbf(coordx,0.)
    zdev_x=0.5*(tempValp-tempValm)/delta
    
    coordy=coordy+delta
    tempValp=rbf(0.,coordy)
    coordy=coordy - 2.*delta
    tempValm=rbf(0.,coordy)
    zdev_y=0.5*(tempValp-tempValm)/delta
    
    return zdev_x, zdev_y

################################################################################
##################      BEST STARTING MEP CALCULATION        ###################
################################################################################  

def bsmep(alat, rbf, n=900, fac=[1.5, 1.5], delta_theta=1):   
    facx=fac[0]
    facy=fac[1]
    alat_x = alat[0]
    alat_y = alat[1]   
    xa =-facx*alat_x
    ya =-facy*alat_y
    xb = facx*alat_x
    yb = facy*alat_y    
    
    g = np.linspace(0,1,n)
    x = (xb-xa)*g+xa
    #delta = (xb - xa) / n
    delta = 0.001 
    
    delta_theta *= np.pi/180
    rep = int(np.pi/delta_theta)    
    data_ss = []
    data_th = []   
    theta = 0
    i = 0
    
    ss_x, ss_y = shear_strength_xy(alat, rbf, bsmep=True)
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
            ss = shear_strength_bsmep(x, y, m, rbf, delta)
            data_th.append(theta)
            data_ss.append(float(ss))
            i +=1
    
    index = data_ss.index(np.amin(np.abs(data_ss)))
    theta = data_th[index]
    ss = data_ss[index]
    
    if theta == np.pi/2. :
        x = np.zeros(n)
        y = (yb-ya)*g+xa
    elif theta == 0. :
        x = (xb-xa)*g+xa
        y = np.zeros(n)
    else:
        x = (xb-xa)*g+xa
        y = np.tan(theta)*x.copy()
    
    print(data_ss)
    # Tiny perturbation of the string
    #x += np.random.rand(len(x))*alat_x/100
    #y += np.random.rand(len(y))*alat_y/100
    
    return x, y, ss, theta

#########################

def shear_strength_bsmep(x, y, m, rbf, delta=0.001): 
    
    zdev=np.zeros(len(x))
#    zdev = derivative_bsmep(rbf, m, x, y, delta)
    
    for i in range(len(x)):
        coordx=x[i]
        coordy=y[i]
        zdev[i]=derivative_bsmep(coordx, coordy, m, rbf, delta)
    
    #Shear strength in GPa
    ss = np.amax(np.abs(zdev))*10.0
    
    return ss

#########################

def derivative_bsmep(coordx, coordy, m, rbf, delta):
    coordx=coordx+delta
    coordy=coordy+m*delta
    tempValp=rbf(coordx,coordy)
    
    coordx=coordx - 2.*delta
    coordy=coordy - 2.*delta*m
    tempValm=rbf(coordx,coordy)
    zdev=0.5*(tempValp-tempValm)/delta
    
    return zdev
