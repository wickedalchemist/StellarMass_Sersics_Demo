#!/astro/data/siesta1/soft/anaconda_2.2.0/bin/python

import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
import ezgal
import glob, re, string, os, pickle



def AG2005_sersics_2dev(clust, r_term, wrapper):
    '''
    Input: Limiting radius for integrating sersic function, stellar population model from ezgal
    Returns: magnitude within limiting radius, equivalent stellar mass, and light to mass ration from ezgal models
    '''

    #LOAD UP AG2005 TABLE 4 VALUES PLUS REDSHIFT 
    dt = np.dtype({'names':['Name','Mtot', 'Mtoterr', 'Mtoterr2', 'Re', 'n', 'z', 'm500', 'm500err', 'M1', 're1', 'M2', 're2', 'mtot', 'mtot_low', 'mtot_high'], 'formats':['S99', np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float]})
    AG = np.loadtxt('/astro/data/siesta1/PROJECTS/ICL/ICL_REDUX/AG2005_Table3n4_complete.txt', dtype=dt)

    #IDENTIFY DESIRED CLUSTER FROM TEXT FILE
    this = np.where(AG['Name']==clust)
    #WE'RE USING THE 2 COMPONENT FITS, SO WE NEED THE MAGNITUDES, EFFECTIVE RADII
    # AND SERSIC INDEX VALUES FOR BOTH COMPONENTS
    Ms = np.array([AG['M1'][this][0],AG['M2'][this][0]],dtype=float)
    Res = np.array([AG['re1'][this][0], AG['re2'][this][0]], dtype=float)
    ns = np.array([4.,4.], dtype=float)
    z = AG['z'][this]

    #ANALYTICAL PERSCRIPTION FOR b_n, THE NORMALIZING FACTOR IN SERSIC PROFILE
    bn = 1.9992*ns-0.3271
    const = np.exp(bn)/(bn**(2*ns))
    #GAMMA FUNCTION FOR EACH SERSIC INDEX, FROM R=O-INF
    gammac = [np.float(mp.gamma(2*n)) for n in ns]

    #USING THE INTEGRATED MAGNITUDES OF EACH COMPONENT WE CAN FIND THE TOTAL INTEGRATED MAGNITUDE
    #WHICH WE NEED TO INTEGRATE THE SERSICS TO A DESIRED RADIUS, AND NOT R=INF
    mue = Ms + 5 * np.log10(Res) + 2.5 * np.log10(2*np.pi*n*const*gammac)
    Fs = 10**(-0.4*(Ms[0]-Ms[1])) #f0 = x*f1
    Minf = Ms[1] + 2.5*np.log10(1/(1+Fs))

    #DESIRED TERMINAL RADIUS
    xterm = bn*(r_term/Res)**(1/ns)
    #INCOMPLETE GAMMA FUNCTION FROM R=0-XTERM
    gammainc = [np.float(mp.gammainc(2*ns[i], a=0, b=xterm[i])) for i in range(len(ns))]
    M = mue - 5 * np.log10(Res) - 2.5 * np.log10(2*np.pi*ns*const*gammainc)
    Fs = 10**(-0.4*(M[0]-M[1])) #f0 = x*f1
    #CAN'T JUST ADD MAGNITUDES, FIND THE SUM OF THE FLUX AND THEN FIND FINAL MAGNITUDE
    Mterm = M[1] + 2.5*np.log10((1+Fs)**(-1))

    #CONVERT TO STELLAR MASS
    if clust=='A0122':
        dm = 38.585890902391228
    else:
        dm = np.squeeze(wrapper.get_distance_moduli(zs=[z], nfilters=1))
    # APPARENT MAGNITUDE OF THE SUN AT THE REDSHIFT OF THE CLUSTER
    msun=np.squeeze(wrapper.get_apparent_mags(2.5, filters='sloan_i', zs=z)) #sloan_i        
    #SOLAR MASS TO LIGHT RATIO FROM EZGAL STELLAR POPULATION MODELS
    mass = np.squeeze(wrapper.get_masses(2.5, zs=z))
    #TOTAL STELLAR MASS WITHIN R=0-INF
    sm_tot = 10**(-0.4*(Minf+dm+0.52-msun))*mass#/1e11
    #TOTAL STELLAR MASS WITHIN DESIRED RADIUS
    sm = 10**(-0.4*(Mterm+dm+0.52-msun))*mass#/1e11)
    
    #IN ADDITION TO RETURNING TERMINAL RADIUS, MAGNITUDE WITHIN DESIRED RADIUS AND STELLAR MASS
    #ALSO RETURN CLUSTER IDENTIFYING CHARACTERISTICS (E.G. MASS, REDSHIFT)
    return r_term, Mterm, sm, AG['z'][this], AG['m500'][this], AG['m500err'][this], mass
      
    
