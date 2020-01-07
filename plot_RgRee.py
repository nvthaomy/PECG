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

plotTitle = 'Rg vs DOP'
ylabel = '$R_g (nm)$' #'$R_g (\sigma)$'
xlabel = 'DOP'
x = [12,24,60,90]
Dirs = [
'xp0.1_N12_f0_V157_LJPME_298K',  
'xp0.1_N24_f0_V157_LJPME_298K',  
'xp0.1_N60_f0_V512_LJPME_298K',  
'xp0.1_N90_f0_V588_LJPME_298K']

dataFName = 'AllStats.dat'
keywords = 'Rg'
avgCol = 0
errCol = 2
lengthScale = 1 #giving lengthScale in Angstrom and Rg calculated by mdtraj will result in Rg in nm
##################
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
    
RgAvgs = np.array(RgAvgs) * lengthScale
RgErrs = np.array(RgErrs) * lengthScale

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.errorbar(x, RgAvgs, yerr=RgErrs, elinewidth=0.75, marker='o', ls=':', lw=1, ms=4, capsize=2)
plt.ylabel(ylabel)
plt.xlabel(xlabel)
title = plotTitle
plt.title(title,loc='center')
plt.savefig('_'.join(title.split())+'.png',dpi=500,bbox_inches='tight')
plt.show()
    
