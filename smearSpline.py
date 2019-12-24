import numpy as np
#from scipy.fft import *
from numpy.fft import *
import spline
import matplotlib, scipy
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


PotName = 'PSplineA_A'
FFfile = 'xp0.1_N12_f0_V157_LJPME_298K_NVT_Spline_ff.dat'
cut = 8.5
#number data points
N = 10000 
#smearing length scale for particle i and j, gamma = 1/(2*np.pi* a**2)**(3./2.) * np.exp(-rs**2/(2* a**2))
ai = 1.5
aj = 1.5
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
    us = B * np.exp(-rs**2/(4*a**2))/(4*np.pi*a**2)**(3./2.)
    return us

# get knots
print('Extracting spline knots from forcefield file')
rs = np.linspace(0,cut,N)
ks = range(0,N)
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
        print ('params')
        print(params)
        for i, param in enumerate(params):
            print('param')
            print(param)
            dict_val = ast.literal_eval('{'+param)
            ParamDict.update(dict_val)
knots = ParamDict['Knots']
print ('Knots for {}:'.format(PotName))
print (knots)

#getting values of potential
print('\nGetting values for pair potential')
us = []
s = spline.Spline(cut, knots)
for r in rs:
    us.append(s.Val(r))

#get smeared density
gammasi = getSmearedDensity(rs, ai)
gammasj = getSmearedDensity(rs, aj)

# Fourier transform gammas and us
print('\nFFT smearing functions and pair potential')
#w = np.blackman(N)
us_fft = fft(us)
gammasi_fft = fft(gammasi)
gammasj_fft = fft(gammasj)
usSmeared_fft = gammasi_fft * gammasj_fft * us_fft
usSmeared = ifft(usSmeared_fft)
#to check FFT
us_ifft = ifft(us_fft)
B = 1. # np.max(usSmeared)
aGauss = np.sqrt(np.mean(ai**2+aj**2))
usGauss = getUGauss(B,rs,aGauss)
#normalized so that maximum potential of Gaussian matches smeared spline
usGauss = usGauss/np.max(usGauss)*np.max(usSmeared)

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
ax.plot(ks, usSmeared_fft, marker = None, ls='-', lw=1, label='FFT$u_{smeared}$')
ax.legend(loc='best',prop={'size': 5})
plt.ylabel('Density')
plt.xlabel('k')
plt.savefig('fftPotential.png',dpi=500,bbox_inches='tight')

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(rs, us, marker = None, ls='-', lw=1, label='no smearing')
ax.plot(rs, usSmeared, marker = None, ls='-', lw=1, label='smeared ai=%3.2f aj=%3.2f'%(ai,aj))
ax.plot(rs, usGauss, marker = None, ls='-', lw=1, label='gaussian')
ax.legend(loc='best',prop={'size': 5})
plt.ylabel('u')
plt.xlabel('r')
plt.savefig('smearedSpline.png',dpi=500,bbox_inches='tight')
plt.show()
print('Gaussian with a below the particle dimension will not match the smeared spline')
