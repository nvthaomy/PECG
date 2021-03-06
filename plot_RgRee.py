""" overlay density profiles from multiple AA MDs of different Uext
need to generate data file with python ~/bin/PEMD/analysis/1d-histogram.py trajectory.dcd topology -axis  -atoms 
 """

from matplotlib import cm, ticker
import numpy as np
import os, sys, glob
import matplotlib
import re
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
################################

plotTitle = 'Ree vs DOP xp0.1 f0 LJPME 298K'
ylabel = '$R_{ee}$ (nm)'
xlabel = 'DOP'
x = [12,24,48,60,90]
Dirs = [
'../xp0.1_N12_f0_V157_LJPME_298K/f0/w0.13/NPT/run2',
'../xp0.1_N24_f0_V157_LJPME_298K/f0/w0.14',
'../xp0.1_N48_f0_V512_LJPME_298K/f0/w0.13',
'../xp0.1_N60_f0_V512_LJPME_298K/f0/w0.14/run1',
'../xp0.1_N90_f0_V588_LJPME_298K/f0/w0.14/run1']

#x = [0,0.17, 0.25, 0.33, 0.5, 0.67, 0.75, 0.83, 0.92, 1.]
#Dirs = [
#'../xp0.1_N12_f0_V157_LJPME_298K/f0/w0.13/NPT/run2',
#'../xp0.1_N12_f0.17_V157_LJPME_298K/f0.17/w0.13/NPT/run2',
#'../xp0.1_N12_f0.25_V157_LJPME_298K/f0.25/w0.13/NPT',
#'../xp0.1_N12_f0.33_V157_LJPME_298K/f0.33/w0.13/run1',
#'../xp0.1_N12_f0.5_V157_LJPME_298K/f0.5/w0.13/NPT/run2',
#'../xp0.1_N12_f0.67_V157_LJPME_298K/f0.67/w0.13/NPT',
#'../xp0.1_N12_f0.75_V157_LJPME_298K/f0.75/w0.13/NPT',
#'../xp0.1_N12_f0.83_V157_LJPME_298K/f0.83/w0.13/NPT',
#'../xp0.1_N12_f0.92_V157_LJPME_298K/f0.92/w0.13/NPT',
#'../xp0.1_N12_f1_V157_LJPME_298K/f1/w0.13/NPT/run2']

dataFName = 'AllStats.dat'
keywords = 'Ree'
avgCol = 0
errCol = 2
fitRg = True
##################
def R(N,b,v):
    return b*N**v

def curvefit(func, x, y):
    from scipy.optimize import curve_fit
    params, params_covariance = curve_fit(func, x, y)
    params_errors = np.sqrt(np.diag(params_covariance))
    return params,params_errors

#check
if len(x) != len(Dirs):
    Exception('Mismatch in sizes of x and Dirs')
RgAvgs = []
RgErrs = []
cwd = os.getcwd()
for i, dir in enumerate(Dirs): 
    num = []
    f = open(os.path.join(dir+'/',dataFName), 'r')
    lines = f.readlines()
    for line in lines:
        if keywords in line:
            vals = line.split()
            for val in vals:
                try:
                    val = float(val)
                    num.append(val)
                except:
                    pass
            RgAvgs.append(num[avgCol])
            RgErrs.append(num[errCol])
if fitRg:    
    [b,v], [bErr,vErr] = curvefit(R, x, RgAvgs)
    str = "$%s = %5.3f \pm %5.3f N^{%5.3f\pm%5.3f}$"%(keywords,b,bErr,v,vErr)
else:
    str = ""
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.errorbar(x, RgAvgs, yerr=RgErrs, elinewidth=0.75, marker='o', ls=':', lw=1, ms=4, capsize=2)
plt.ylabel(ylabel)
plt.xlabel(xlabel)
plt.text(0.75*np.mean(x),np.mean(RgAvgs),str)
title = plotTitle
plt.title(title,loc='center')
plt.savefig('_'.join(title.split())+'.png',dpi=500,bbox_inches='tight')
plt.show()
    
