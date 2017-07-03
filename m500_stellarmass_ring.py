#!/astro/data/siesta1/soft/anaconda_2.2.0/bin/python

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import mpmath as mp
import ezgal
import scipy.odr as odr
import glob, re, string, os, pickle, glob

import good_clusters
import cluster_class
import perl_sb
import integrate_sersics

#######################
# ODR FITTING FUNCTION(S) DEFINITIONS
######################
def lin_func(p,x):
    m, b = p
    return m*x+b

def fixed_lin_func(intercept,x):
    #print "this is b: ", intercept
    return fixed_slope*x+intercept

def ODR(x,y,x_err=False,y_err=False):
    #A LA http://stackoverflow.com/questions/22670057/linear-fitting-in-python-with-uncertainty-in-both-x-and-y-coordinates
    
    #Create a model for fitting
    lin_model = odr.Model(lin_func)
    #Set up RealData object
    data = odr.RealData(x,y,sx=x_err, sy=y_err)
    #Set up ODR with model and data
    tmp = odr.ODR(data, lin_model, beta0=[0.,1.])
    #Run regression
    out = tmp.run()
    #out.pprint()
    return out

####################
# FRACTIONAL UNCERTAINTY ON MEASURED MAGNITUDES
####################

def frac_magerr(mtot, mtot_low, mtot_high, smterm, sol_mass):
    '''
    Given an uncertainty on a total magnitude (r<inf), find the fractional error
    And return uncertainty on a magnitude within a specified radius with the same fractional errror
    '''
    frac_low = 1 - 10**(-0.4*(-1*mtot_low))
    frac_high = 1 - 10**(-0.4*(mtot_high))

    dlow = -2.5*np.log10(smterm/(frac_low*smterm+smterm))
    dhigh = -2.5*np.log10(smterm/(frac_high*smterm+smterm))

    dMlow = sol_mass*0.4*np.log(10)*smterm*dlow
    dMhigh = sol_mass*0.4*np.log(10)*smterm*dhigh

    sm_err = 0.5 * (smterm+dMhigh-(smterm-dMlow))

    return sm_err

####################
# SET UP FIGURE
###################

fig = plt.figure(1, figsize=(8,8))
ax = fig.add_subplot(111)


###################
# CHOICE OF IMF
###################
wrapper=ezgal.wrapper(glob.glob( '/astro/data/siesta1/PROJECTS/ICL/ICL_REDUX/SPS_models/*bc03_ssp_z_0.02_chab.model' ))

####################
# CHOOSE RING RADII
####################

term2 = 100
term1 = 10

#####################
# AG2005 SAMPLE
#####################
dt = np.dtype({'names':['Name','Mtot', 'Mtoterr', 'Mtoterr2', 'Re', 'n', 'z', 'm500', 'm500err', 'M1', 're1', 'M2', 're2', 'mtot', 'mtot_low', 'mtot_high'], 'formats':['S99', np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float]})
AG = np.loadtxt('/astro/data/siesta1/PROJECTS/ICL/ICL_REDUX/AG2005_Table3n4_complete.txt', dtype=dt)
#ALL M500 MASSES ARE SYSTEMATICALLY LOWER THAN OTHER SAMPLES BY ~15%
AG['m500']=AG['m500']*1.15

AG_SM = []
AG_SMerr = []
for i, clust in enumerate(AG['Name']):
    #USE INTEGRATE_SERSIC FUNCITON TO GET STELLAR MASS WITHIN AN ANULUS WITH SPECIFIED INNER,OUTER RADII
    r_term, Mterm, smterm2, z, m500, m500err, sol_mass = integrate_sersics.AG2005_sersics_2dev(clust, term2, wrapper)
    r_term, Mterm, smterm1, z, m500, m500err, sol_mass = integrate_sersics.AG2005_sersics_2dev(clust, term1, wrapper)

    #GET UNCERTAINTY ON STELLAR MASS USING FRAC_MAGERR, DEFINED ABOVE
    dM1 = frac_magerr(AG['mtot'][i], AG['mtot_low'][i], AG['mtot_high'][i], smterm1, sol_mass)
    dM2 = frac_magerr(AG['mtot'][i], AG['mtot_low'][i], AG['mtot_high'][i], smterm2, sol_mass)

    if term1 != 0:
        AG_SM.append(smterm2-smterm1)
        AG_SMerr.append(np.sqrt(dM1**2+dM2**2))
    else:
        AG_SM.append(smterm2)
        AG_SMerr.append(dM2)

#PLOT STELLAR MASS VS CLUSTER MASS (M500), COLOR CODED BY REDSHIFT
ax.errorbar(AG['m500']*1e14, AG_SM, yerr = [AG_SMerr, AG_SMerr], xerr=[AG['m500err']*1e14, AG['m500err']*1e14], color='black', marker = 'None',ls='None', alpha=0.75)
cs = cm.coolwarm(AG['z']/(1.75))
ax.scatter(AG['m500']*1e14, AG_SM, color=cs, s=200, marker='*')

########################
# GET STELLAR MASSES & M500 FROM OTHER SAMPLES PREVIOUSLY DEFINED AND STORED IN A PICKLE
#########################

pkl_name = 'smring_'+str(term1)+'-'+str(term2)+'.pkl'
pkl = pickle.load(open(pkl_name, 'rb'))
m500 = pkl['m500']
m500_err = pkl['m500_err']
sm = pkl['sm']
sm_err = pkl['sm_err']
kind = pkl['kind']
clusters = pkl['clusters']
zs = pkl['zs']

############################
# DEFINE SUB-SAMPLES AND FIT WITH ODR
#############################

#DEFINE SUB-SAMPLES
low_z = np.where( m500[kind=='AG']>0 )
mid_z = np.where((clusters != 'CL1226') & (kind=='PII'))
high_z = np.where( (kind=='PERL') )

#FIT THE MID-REDHSIFT SAMPLE RELATION LOG(STELLAR MASS) = ALPHA*LOG(M500)+BETA WITH ODR
x = np.linspace(np.log10(1e13), np.log10(3e15))
logx = np.log10(m500[mid_z])
xerr = np.array([0.5*(np.log10(m500[mid_z][i]+m500_err[mid_z][i][1])-np.log10(m500[mid_z][i]-m500_err[mid_z][i][0])/2e14) for i in range(len(m500[mid_z]))])
logy = np.log10(sm[mid_z])
yerr = np.log10(sm_err[mid_z])
pii = ODR(logx, logy, y_err=yerr, x_err=xerr)
fixed_slope = pii.beta[0]
#PLOT MID-Z SAMPLE
ax.errorbar(m500[mid_z], sm[mid_z], yerr = [sm_err[mid_z], sm_err[mid_z]], xerr=10**xerr, color='navy', marker = 'o', ms=8, lw=2, ls='None', alpha=0.75)

#FIT HIGH-REDSHIFT SAMPLE WITH SAME RELAITON, ASSUMING SAME SLOPE AS FOR MID-REDSHIFT SAMPLE (AS FIT ABOVE)
logx = np.log10(m500[high_z])
xerr = np.array([0.5*(np.log10(m500[high_z][i]+m500_err[high_z][i][1])-np.log10(m500[high_z][i]-m500_err[high_z][i][0])/2e14) for i in range(len(m500[high_z]))])
logy = np.log10(sm[high_z])
yerr = np.log10(sm_err[high_z])
#Create a model for fitting
fixed_lin_model = odr.Model(fixed_lin_func)
#Set up RealData object
data = odr.RealData(logx,logy,sx=xerr, sy=yerr)
#Set up ODR with model and data
tmp = odr.ODR(data, fixed_lin_model, beta0=[5.])
#Run regression
hz = tmp.run()
#PLOT HIGH-Z SAMPLE
ax.errorbar(m500[high_z], sm[high_z], yerr = sm_err[high_z], xerr=10**xerr, color='red', marker = 'D', ms=8, lw=2, ls='None', alpha=0.75)

#PLOT BEST-FIT RELATIONS
x = np.linspace(np.log10(1e13), np.log10(3e15))
ax.plot(10**x, 10**(fixed_slope*x+pii.beta[1]), lw=2, ls='-', color='navy')
ax.plot(10**x, 10**(fixed_slope*x+hz.beta), lw=2, ls='--', color='tomato')

##########
# make points off plot for legend
#########
ax.plot(-99, -99, color='k', marker='*', ms=15, lw=2, ls='None', label='GZZ13' )
ax.plot(-99, -99, color='k', marker='^', ms=10, lw=2, ls='None', label='Paper II' )
ax.plot(-99, -99, color='k', marker='D', ms=10, lw=2, ls='None', label='high-z' )

plt.semilogy()
plt.semilogx()
plt.xlim(1.4e13, 1.1e15)
plt.ylim(9e10, 2.5e12)

#axis line width
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)

#ticks
#ax.minorticks_on()
ax.tick_params('both', length=15, width=2, which='major', labelsize=14)
ax.tick_params('both', length=8, width=1, which='minor')

#labels
plt.xlabel(r'M$_{500}$ [M$_\odot$]', fontsize=18, fontweight='bold')
plt.ylabel(r'M$_\bigstar$ [M$_\odot$]', fontsize=18, fontweight='bold')
plt.title(str(int(term1))+' < r [kpc] < '+str(int(term2)), fontsize=18, fontweight='bold')
plt.legend(loc = 'upper left', numpoints = 1)

plt.savefig('sm'+str(int(term2))+'-'+str(int(term1))+'_m500.png')

plt.show()
plt.clf()



