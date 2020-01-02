import numpy as np
#from scipy.fft import *
from numpy.fft import *
import spline
import matplotlib, scipy, cmath
from scipy import integrate
import re, os, ast
showPlots = True
try:
  os.environ["DISPLAY"] #Detects if display is available
except KeyError:
  showPlots = False
  matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
import matplotlib.pyplot as plt
#plt.style.use('seaborn-dark')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
colors = ['#6495ED','r','#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#1E90FF']
#################
"""Smear spline pair potential"""

PotName = 'PSplineA-_Na+'
FFfile = 'xp0.1_N12_f1_V157_LJPME_298K_NVT_Spline_ff.dat'
OutFF = 'smeared_'+FFfile
cut = 8.5
#number data points
N = 2000 
#number of knots for smeared spline
Nknots = 10
#smearing length scale for particle i and j, gamma = 1/(2*np.pi* a**2)**(3./2.) * np.exp(-rs**2/(2* a**2))
ai = 1. 
aj = 1.
#smearing length scale of Gaussian, u = 1/(4*np.pi* aij**2)**(3./2.) * np.exp(-rs**2/(4 * aij**2)), aij**2 = (ai**2+aj**2)/2
aiG = 1.
ajG = 1.
#----------
def getSmearedDensity(rs, a):
    rs = np.array(rs)
    a = float(a)
    gammas = np.exp(-rs**2/(2* a**2))
    norm = np.sum(gammas)
    gammas /= norm    
    return gammas

def getUGauss(B,rs,a):
    rs = np.array(rs)
    a = float(a)
    B = float(B)
    us = B * np.exp(-rs**2/(2*a**2))/(2*np.pi*a**2)**(3./2.) 
#    us = B * np.exp(-rs**2/(4*a**2))/(4*np.pi*a**2)**(3./2.)
    return us
def getExVol(rs, us):
    us = np.array(us)
    rs = np.array(rs)
    rs2 = rs**2
    mayerf =  np.exp(-us) - 1
    y =  4 * np.pi * rs2 * (-mayerf)
    v = integrate.simps(y, rs)
    return v
def extrap(usSmeared, method = 'MaxSlope'):
    """method: MaxSlope: do linear extrapolation in the region of r = 0 to r where slope is maximum
               ZeroSlope: do linear extrapolation with slope = 0 in the region of r = 0 to r of the first peak"""
    usSmeared_ext = usSmeared.copy()
    usSmeared_max = np.max(usSmeared)
    usSmeared_max_id = np.where(abs(usSmeared-usSmeared_max)<1e-6)[0][0]
    if method == 'MaxSlope':
        #linear extrapolate from r = 0 to r where du/dr is maximum
        du_dr = (usSmeared[:-1] - usSmeared[1:])/(rs[:-1]-rs[1:])
        #find r where abs(du/dr) is maximum, only find max slope beyond rmax
        du_dr_max = np.max(np.abs(du_dr[usSmeared_max_id:]))
        du_dr_max_id = np.where(abs(abs(du_dr[usSmeared_max_id:])-du_dr_max)<1e-10)[0][0]
        du_dr_max_id += usSmeared_max_id
        #extrapolate inner region of smeared spline
        du_dr_max = du_dr[du_dr_max_id]
        usSmeared_ext[:du_dr_max_id] = -du_dr_max * (rs[du_dr_max_id]-rs[:du_dr_max_id]) + usSmeared[du_dr_max_id]
    elif method == 'ZeroSlope':
        usSmeared_ext[:usSmeared_max_id] = [usSmeared_max]*len(range(0,usSmeared_max_id))
    return usSmeared_ext
   
# get knots
print('\n\t===Extracting spline knots from forcefield file==')
rs = np.linspace(0,cut,N/2)
knots = []
f = open(FFfile, 'r')
str = f.read()
f.close()
ParamDict = {} #dictionary of paramters
potential_str = [val for val in str.split('>>> POTENTIAL ') if len(val) > 0]
#reading ff file and extracting knot values
for potential in potential_str:
    s = potential.split('{') #split potential name from params
    if  PotName in s[0]:
        s = s[-1] #only take the parameters
        params = [s] #[' '.join(val.rsplit()) for val in re.split('}|,',s) if len(val)>0] 
        for i, param in enumerate(params):
            dict_val = ast.literal_eval('{'+param)
            ParamDict.update(dict_val)
knots = ParamDict['Knots']
print ('Knots for {}:'.format(PotName))
print (knots)

#getting values of potential
print('\n\t===Getting values for pair potential===')
us = []
s = spline.Spline(cut, knots)
for r in rs:
    us.append(s.Val(r))
# Making symmetrical data around r = 0
rs_rev = -1*rs[::-1][:-1]
rs = np.array(rs_rev.tolist() + rs.tolist())
ks = range(0,len(rs))
if not isinstance(us, list):
    us = us.tolist()
us =  np.array(us[::-1][:-1] + us)

#get smeared density
gammasi = getSmearedDensity(rs, ai)
gammasj = getSmearedDensity(rs, aj)

# Fourier transform gammas and us
print('\n\t===FFT smearing functions and pair potential===')
w = np.ones(len(rs))
#w = np.blackman(N)
us_fft = fft(us*w)
gammasi_fft = fft(gammasi*w)
gammasj_fft = fft(gammasj*w)
usSmeared_fft = gammasi_fft * gammasj_fft * us_fft
usSmeared = ifft(usSmeared_fft)
#only take the real part (imaginary is close to 0)
usSmeared = usSmeared.real

#to check FFT
us_ifft = ifft(us_fft)


#extrapolate spline to r=0
#print('\n\t===Extrapolating spline to r = 0===')
#usSmeared_ext = extrap(usSmeared, method = '')
usSmeared_ext = usSmeared

# get new spline knots
print('\n\t===Getting new spline knots===')
knots1 = Nknots * [0.]
r0_id = np.where(abs(rs-0) < 1e-5)[0][0]
s1 = spline.Spline(cut, knots1)
s1.fitCoeff(rs[r0_id:], usSmeared_ext[r0_id:])
knots1 = s1.knots.tolist()
#print ("New knots: \n{}".format(knots1))

#write ff data
for i, val in enumerate(potential_str):
    if PotName in val: #find index of the old param string and remove it
        potential_str[i] = ''
SplineStr = []
s = PotName
s += "\n{"
s += "'Knots' : {} ".format(knots1)
s += "}\n"
potential_str[i] = s
str = [val for val in potential_str if len(val) > 0]
str = '>>> POTENTIAL '.join(str)
str = '>>> POTENTIAL ' + str
ff = open(OutFF,'w')
ff.write(str)

#get Gaussian potential
B = 1. # np.max(usSmeared)
aGauss = np.sqrt((aiG**2+ajG**2))
usGauss = getUGauss(B,rs,aGauss)
#normalized so that maximum potential of Gaussian matches smeared spline
usGauss = usGauss/np.max(usGauss)*np.max(usSmeared_ext)

#plot
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(rs,gammasi, marker = None, ls='-', lw=1, label='$\Gamma_i$')
ax.plot(rs,gammasj, marker = None, ls=':', lw=1, label='$\Gamma_j$')
ax.legend(loc='best',prop={'size': 5})
plt.ylabel('Density')
plt.xlabel('r')
plt.savefig('smearingFunctions.png',dpi=500,bbox_inches='tight')

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(ks, gammasi_fft, marker = None, ls='-', lw=1, label='FFT $\Gamma_i$')
ax.plot(ks, gammasj_fft, marker = None, ls=':', lw=1, label='FFT $\Gamma_i$')
ax.legend(loc='best',prop={'size': 5})
plt.ylabel('Density')
plt.xlabel('k')
plt.savefig('fftSmearingFunc.png',dpi=500,bbox_inches='tight')

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(ks, us_fft, marker = None, ls='-', lw=1, label='FFT u')
ax.plot(ks, gammasi_fft, marker = None, ls='-', lw=1, label='FFT $\Gamma_i$')
ax.plot(ks, gammasj_fft, marker = None, ls=':', lw=1, label='FFT $\Gamma_i$')
ax.plot(ks, usSmeared_fft, marker = None, ls='-', lw=1, label='FFT $u_{smeared}$')
ax.legend(loc='best',prop={'size': 5})
#plt.ylabel('Density')
plt.xlabel('k')
plt.savefig('fftPotential.png',dpi=500,bbox_inches='tight')


fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(rs, us, marker = None, ls='-', lw=1, label='no smearing')
ax.plot(rs, usSmeared, marker = None, ls='-', lw=1, label='smeared ai=%3.2f aj=%3.2f'%(ai,aj))
ax.plot(rs, usGauss, marker = None, ls='-', lw=1, label='gaussian aiG=%3.2f ajG=%3.2f'%(aiG,ajG))
ax.legend(loc='best',prop={'size': 5})
plt.ylabel('u')
plt.xlabel('r')
plt.xlim(0)
plt.savefig('smearedSpline.png',dpi=500,bbox_inches='tight')
#plt.show()
print('\n===Gaussian with a below the particle dimension will not match the smeared spline===')

# get excluded volume
vGauss = getExVol(rs, usGauss)
vSpline = getExVol(rs, us)
vSmearedSpline = getExVol(rs, usSmeared_ext)
print(vSpline)
vs = np.array([vGauss, vSpline, vSmearedSpline])
radii = (vs * 3. /(4*np.pi))**(1./3.)
print('\n\tRadii estimated from excluded volume for Gauss Spline SmearedSpline:')
print('%3.4f %3.4f %3.4f'%(radii[0],radii[1],radii[2]))
