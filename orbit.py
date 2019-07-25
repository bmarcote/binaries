#!/usr/bin/python

##########
# NOTE: This is script from Javi that I need to look at. It is not part of the
# project. And I have a different version.
#
#


import numpy as np
import sys
from pylab import *
import matplotlib
from matplotlib import*
import matplotlib.pyplot as plt
from cmath import *
from scipy.optimize import fsolve
from scipy.special import gamma as g_func
from scipy import integrate
from scipy import savetxt, transpose
from matplotlib.patches import Ellipse
from matplotlib.ticker import FormatStrFormatter,MaxNLocator,MultipleLocator



twocolumn = False

if twocolumn == True:
    fig_width_pt = 255.76535
    fontzise     = 8
else:
    fig_width_pt = 426.79134
    fontzise     = 10
    
    
#fig_width_pt = 255.76535  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt = 426.79134  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)+1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width/golden_mean      # height in inches
fig_height = fig_height
#fig_height = fig_width*np.sqrt(1.6)
fig_size =  [fig_width,fig_height]

params = {'backend': 'eps',
		'axes.labelsize': fontzise,
		'text.fontsize': fontzise,
		'legend.fontsize': fontzise,
		'xtick.labelsize': fontzise,
		'ytick.labelsize': fontzise,
		'axes.linewidth': 0.5,
		'linewidth': 0.5,
		'text.usetex': True,
		'ps.usedistiller': False,
		'figure.figsize': fig_size,
		'font.family': 'Times New Roman',
		'font.fantasy': ['Comic Sans MS'],
		'font.serif': ['Bitstream Vera Serif']}
pylab.rcParams.update(params)




pi = np.pi
radgra = 180/pi
grarad = pi/180
G = 6.6732e-8
AU = 1.49597892e13
MS = 1.989e33
RScm = 6.955e10
RS = RScm/AU


parametros = open('orb_param.par', 'r')
lines = parametros.readlines()
name      = lines[3][33:42]                # Nombre del objeto           
M1        = float(lines[4][33:42])         # Masa de la primaria         
R1        = float(lines[5][33:42])         # Radio de la primaria
M2        = float(lines[6][33:42])         # Masa del objeto compacto    
P         = float(lines[7][33:42])         # Periodo orbital             
e         = float(lines[8][33:42])         # Excentricitat               
phase_per = float(lines[9][33:42])         # Fase del periastro          
wp        = float(lines[10][33:42])         # Argumento del periastro     
lan       = float(lines[11][33:42])        # Longitud del nodo ascendente
inc       = float(lines[12][33:42])        # Inclinacion de la orbita    
n         = float(lines[13][33:42])        # Numero de puntos            
phas_ini  = float(lines[14][33:42])        # Fase inicial                
phas_fin  = float(lines[15][33:42])        # Fase final                  
phas_infc = float(lines[16][33:42])        # Fase inferior conjunction                
phas_supc = float(lines[17][33:42])        # Fase superior conjunction   

parametros.close()

M1=M1*MS
M2=M2*MS
R1=R1*RScm

Psec = P*86400.0

mu = G*(M1+M2)
a = (mu*Psec**2.0/(4.0*pi*pi))**(1.0/3.0)


def position(phase, inc=1.0):
    """
    position(MJD):
        Return the 'x' (along the LOS with + being towards us) and 'y' (in the
            plane of the sky with + being away from the line of nodes and -
            being in the direction of the line of nodes) positions of the
            pulsar with respect to the center of mass in units of lt-sec.
            (Note:  This places the observer at (+inf,0.0) and the line of nodes
            extending towards (0.0,-inf) with the pulsar orbiting (0.0,0.0)
            counterclockwise).  'inc' is the inclination of the orbit in degrees.
            MJD can be an array.  The return value is (xs, ys).
    """
    ma, ea, ta = calc_anoms(phase)
    ws = wp*grarad
    orb_phs = ta + ws
    sini = np.sin(inc*grarad)
    r = a*(1.0-e*e)/(1.0+e*np.cos(ta))
    x = r*np.sin(orb_phs)*np.sin(inc*grarad)
    y = -r*np.cos(orb_phs)
    z =  r*np.sin(orb_phs)*np.cos(inc*grarad)
    vx = -2*pi*a*a*np.sin(ea)/(r*Psec)
    vy =  2*pi*a*a*np.sqrt(1-e*e)*np.cos(ea)/(r*Psec)
    vrel = np.sqrt(vx*vx+vy*vy)
    return x, y, z, vrel, ta

def calc_anoms(phase):
    """
    calc_anoms(phase):
        Return a tuple of the mean, eccentric, and true anomalies (all in radians)
            for every phase.
    """
    mean_anom = mean_anomaly(phase)
    ecc_anom = eccentric_anomaly(mean_anom)
    true_anom = true_anomaly(ecc_anom, e)
    return (mean_anom, ecc_anom, true_anom)

def mean_anomaly(phase):
    return 2*pi*(phase - phase_per)

def eccentric_anomaly(mean_anomaly):
    """
    eccentric_anomaly(mean_anomaly):
        Return the eccentric anomaly in radians, given a set of mean_anomalies
        in radians.
    """
    ma = np.fmod(mean_anomaly, 2*pi)
    ma = np.where(ma < 0.0, ma+2*pi, ma)
    eccentricity = e
    ecc_anom_old = ma
    ecc_anom = ma + eccentricity*np.sin(ecc_anom_old)
    # This is a simple iteration to solve Kepler's Equation
    while (np.maximum.reduce(np.fabs(ecc_anom-ecc_anom_old)) > 5e-15):
        ecc_anom_old = ecc_anom[:]
        ecc_anom = ma + eccentricity*np.sin(ecc_anom_old)
    return ecc_anom

def true_anomaly(ecc_anom, e):
    """
    true_anomaly(ecc_anom, e):
        Return the True Anomaly (in radians) given the Eccentric anomaly
            (ecc_anom in radians) and the eccentricity (e)
    """
    return 2.0*np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(ecc_anom/2.0))

def rot_lan(x,y,z):
    """
    rot_lan(x,y,z):
        Return the coordinates for every point of the orbit rotated by the
        longitude of the ascending node (lan). Definido el origen en la 
    """
    xp = x
    yp =  y*np.cos(lan*grarad) + z*np.sin(lan*grarad)
    zp = -y*np.sin(lan*grarad) + z*np.cos(lan*grarad)
    return xp, yp, zp


# ***  Main program  ***

phase = np.linspace(phas_ini,phas_fin,n+1)    # 1. Generating the orbital phases

phase_conjunctions = np.array([phas_infc,phas_supc])# 1.1 Generating the conjunction orbital phases

phase_periast = np.array([phase_per, phase_per+0.5])# 1.2 Generating periastron and apastron orbital phases

phase_points = np.linspace(0.0,1.0, 11) # 1.3 Generating points to mark the position of each 0.1 orbital phase

x, y, z, vrel, theta = position(phase, inc)          # 2. Computing the orbit
xp, yp, zp, vrelp, thetap = position(phase_points, inc)
xc,yc,zc,vrelc,thetac = position(phase_conjunctions, inc)
xper, yper, zper, vrelper, thetaper = position(phase_periast, inc)

factor = AU

linewidth_obs = 4.0

#print x, y, z, theta

fig = plt.figure()

ax = fig.add_subplot(111, aspect='equal')

# Plotting the star
star = Ellipse(xy=(1.e-30,1.e-30), width = 2*R1/factor, height=2*R1/factor, angle = 0., ec='darkgrey',fc='darkgrey', fill=True,zorder=100,alpha=1.)
ax.add_artist(star)

#plotting the elliptical orbit
ax.plot(z/factor,y/factor,'-k')
# plotting the 10 circular points in 0.1 orbital phase steps
ax.plot(zp/factor,yp/factor,marker='o', mec='k', mfc='w', lw=0.0, ms=5., zorder=10)
# plotting the inferior and superior conjunctions
ax.plot(zc/factor,yc/factor, marker='s', mec='k', mfc='w', color='k', ms=5., linestyle='dotted', lw=0.0, zorder=1000)
# plotting position of periastron and apastron
ax.plot(zper/factor,yper/factor, marker='v', mec='k', mfc='w', color='b',ms=5., linestyle='dotted', lw=0.5, zorder=100)


# Annotating useful information into the plot:

ax.annotate('', xy=(1.1, 0.0), xytext= (1.1,0.25), arrowprops=dict(arrowstyle="-|>", linewidth=1., shrinkA=2., shrinkB=0.,fc = 'k', ec = 'k',mutation_scale=8.), fontsize=8)
ax.annotate('To observer', xy=(1.1, 0.0), xytext= (0.92,-0.1), fontsize=10)
ax.annotate('$\\phi=0$', xy=(0.7, -0.4), xytext= (0.7,-0.4), fontsize=10)

ax.annotate('$\\phi_{ic}=0.84$', xy=(0., -0.9), xytext= (0.,-0.9), fontsize=10)
ax.annotate('$\\phi_{sc}=0.28$', xy=(0., 0.9), xytext= (0.,0.9), fontsize=10)


# Plotting the curved arrow to indicate the orbital motion direction

ax.annotate('', xy=(-0.82, 0.6),  xycoords='data', xytext=(-0.6, 0.85), size=20, arrowprops=dict(arrowstyle="-|>", linewidth=1., shrinkA=2., shrinkB=0.,fc = 'k', ec = 'k',mutation_scale=8., connectionstyle="arc3,rad=0.25"))
            #bbox=dict(boxstyle="round", fc="0.8"),



# Setting the limits and labels:

ax.set_xlim(-1.0,2.2)
ax.set_ylim(-0.999,1.2)
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.5))

ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))

ax.set_xlabel('AU')
ax.set_ylabel('AU')

leg= ax.legend(loc=(0.52,0.62), numpoints=1, markerscale=0.5)
#leg.draw_frame(False)


fig.savefig('orbita.eps', bbox_inches='tight', padinches=0.001, dpi=300)
