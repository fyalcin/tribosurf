import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm

################################################################################
###########################      PLOTTING PDF       ###########################
################################################################################

def plot_pes(alat, xnew, ynew, Enew, element='interface', miller=''):
    alat_x = alat[0]
    alat_y = alat[1]
    fact=1.
    level= 43
    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    anglerot='vertical'
    shrin=1.
    zt1=plt.contourf(xnew, ynew, Enew, level, extent=(-fact*alat_x,fact*alat_x,-fact*alat_y,fact*alat_y), cmap=plt.cm.RdYlBu_r)
    cbar1=plt.colorbar(zt1,ax=ax,orientation=anglerot,shrink=shrin)
    cbar1.set_label(r'$E_{adh} (J/m^2)$', rotation=270, labelpad=20,fontsize=15,family='serif')
    plt.title("PES for " + element + miller, fontsize=18,family='serif')
    ax.quiver(0. , 0., 1., 0.,scale=1.,scale_units='inches',width=0.01,color='white')
    ax.quiver(0. , 0., 0., 1.,scale=1.,scale_units='inches',width=0.01,color='white')
    ax.plot(0.,0.,'w.',ms=7)
    ax.text(0.5,-0.5,'[1 0 1]',rotation='horizontal',color='white', fontsize=14)
    ax.text(-0.5,1.,'[1 2 1]',rotation='vertical',color='white', fontsize=14)
    ax.axis([-fact*alat_x,fact*alat_x,-fact*alat_y,fact*alat_y])
    plt.xlabel(r"distance ($\AA$)",fontsize=12,family='serif')
    plt.ylabel(r"distance ($\AA$)",fontsize=12,family='serif')
    
    for zt1 in zt1.collections:
       zt1.set_edgecolor("face")
       zt1.set_linewidth(0.000000000001)
    
    plt.savefig("PES" + str(element) + str(miller) + ".pdf")

def plot_pep_xy(alat, rbf, element='interface', miller='', npoints=300, fact=1.):
    
    alat_x = alat[0]
    alat_y = alat[1]
    x=np.arange(-1.5*alat_x, 1.5*alat_x, alat_x/npoints)
    y=np.arange(-1.5*alat_y, 1.5*alat_y, alat_y/npoints)
    
    fig1, ( (s1,s2) ) = plt.subplots(nrows=2, ncols=1, figsize=(3, 6), dpi=100, sharex=False, sharey=True)
    fig1.subplots_adjust(left=0.25, bottom=0.12, right=0.99, top=0.95, wspace=0.25, hspace=0.25)
    plt.minorticks_on()
    s1.plot(x, rbf(x,np.zeros(len(x))),'-',color='blue')  
    s1.set_xlim( (-fact*alat_x, fact*alat_x) )
    s1.set_xlabel(r"distance along x ($\AA$)",fontsize=12,family='serif')
    s1.set_ylabel(r'$E_{adh} (J/m^2)$',fontsize=12,family='serif')
    s2.plot(y, rbf(np.zeros(len(y)), y),'-',color='blue')
    s2.set_xlim( (-fact*alat_y, fact*alat_y) )
    s2.set_xlabel(r"distance along y ($\AA$)",fontsize=12,family='serif')
    s2.set_ylabel(r'$E_{adh} (J/m^2)$',fontsize=12,family='serif')
    plt.savefig("PEPxy" + str(element) + str(miller) + ".pdf")
    
def plot_mep(alat, rbf, x, y, xnew, ynew, theta, n=101, element='interface', lattice='', miller='', fac=[0.75, 0.75], fact=1.):   
    facx=fac[0]; facy=fac[1]
    alat_x = alat[0]
    alat_y = alat[1]
    
    xa =-facx*alat_x; ya =-facy*alat_y
    xb = facx*alat_x; yb = facy*alat_y    
    g = np.linspace(0,1,n)

    if theta == 90 or theta == 270: # y direction
        xi = np.zeros(n)
        yi = (yb-ya)*g+ya
    else: # best fing direction or x direction
        xi = (xb-xa)*g+xa
        yi = np.tan(theta)*x.copy()
    
    fig3 = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig3.add_subplot(111)
    ax.set_aspect('equal')
    anglerot='vertical'
    shrin=1.
    level=43
    
    zt1=plt.contourf(xnew, ynew, rbf(xnew, ynew), level, extent=(-fact*alat_x,fact*alat_x,-fact*alat_y,fact*alat_y), cmap=plt.cm.RdYlBu_r)
    cbar1=plt.colorbar(zt1,ax=ax,orientation=anglerot,shrink=shrin)
    cbar1.set_label(r'$E_{adh} (J/m^2)$', rotation=270, labelpad=20,fontsize=15,family='serif')
    
    plt.title("PES/MEP for " + str(element) + str(lattice) + str(miller),fontsize=18,family='serif')
    ax.plot(xi, yi,'.-', c='white', ms=1)
    ax.plot(x, y,'.-', c='yellow', ms=2)
    ax.axis([-fact*alat_x,fact*alat_x,-fact*alat_y,fact*alat_y])
    plt.xlabel(r"distance ($\AA$)",fontsize=12,family='serif')
    plt.ylabel(r"distance ($\AA$)",fontsize=12,family='serif')
    plt.savefig("mep"+ str(element) + str(lattice) + str(miller) +".pdf")


def plot_mep1d(data, element='interface', lattice='', miller=''):
    lxy, dVx, dVy, Vz, Ve, force = data.T
    fig = plt.figure(figsize=(10, 9), dpi=100)
    a1 = fig.add_subplot(221)
    a1.plot(lxy,Vz,'g-',lxy, Ve,'b.',ms=3)
    a2 = fig.add_subplot(222)
    a2.plot(lxy, 10.*np.sqrt(dVx**2 + dVy**2),'g-',lw=0.5)
    a2.plot(lxy, 10.*force,'.-',color='blue',ms=3)
    plt.savefig("mep-1d" + str(element) + str(lattice) + str(miller) +".pdf")
    
################################################################################
###########################      WRITE ON FILE       ###########################
################################################################################

def save_cell_info(out='output.txt', alat=[], scale=0.):
    expansion = str(scale*100)+'%'
    with open(out, 'a') as out:
        out.write("alat_x: {0}\t alat_y: {1}\t lateral expansion: {2}\n".format(str(alat[0]),str(alat[1]), str(expansion)))

def save_ss_xy(out='output.txt', ss_x=None, ss_y=None, element='interface', miller=''):
    with open(out, 'a') as out:
        out.write("")
        out.write("### {0} ({1})\n".format(str(element),str(miller)))
        out.write("Shear strength estimated along x  =  {0}  GPa\n".format(str(ss_x)))
        out.write("Shear strength estimated along y  =  {0}  GPa\n".format(str(ss_y)))

def save_bsmep(out='output.txt', theta=None, data=[]):
    with open(out, 'a') as out:
        out.write("Best starting string for MEP: angle {0} with x axis\n".format(str(theta*180/np.pi))) 
    np.savetxt('BSMEP.txt', data)

def starting_mep(out='output.txt'):
    with open(out, 'a') as out:
        out.write("\nStarting (local) calculation of MEP. This might take a while...\n")

def save_mep(out='output.txt', tol1=1e-7, tol=None, n_nodes=None, nstep=None, data=[]):
    with open(out, 'a') as out:
        out.write("ZTS calculation x with {0} images\n".format(str(n_nodes)))
        if tol > tol1 :
            with open('output_correct.txt', 'a') as out:
                out.write("The calculation failed to converge after {0} iterations tol = {1}\n".format(str(nstep), str(tol)))
        else:
            with open('output_correct.txt', 'a') as out:
                out.write("The calculation terminated after {0} iterations tol = {1}\n".format(str(nstep), str(tol)))
    np.savetxt('MEP.txt', data)
    
def save_ss(out='output.txt', ss_max=None, ss_min=None, data=[]):
    with open(out, 'a') as out:
        out.write("Shear strength estimated along MEP from min/max force/area = {0} / {1}  GPa\n".format(str(ss_min),str(ss_max)))
    np.savetxt('DATA_SS_MEP.txt', data)
