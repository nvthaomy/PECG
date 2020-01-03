""" overlay density profiles from multiple AA MDs of different Uext
need to generate data file with python ~/bin/PEMD/analysis/1d-histogram.py trajectory.dcd topology -axis  -atoms 
 """

from matplotlib import cm, ticker
import numpy as np
import os, sys, glob
import matplotlib
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
################################

legends = ['UextNaCl=0.5', 'UextNaCl=1', 'UextNaCl=2', 'UextNaCl=2.5', 'UextNaCl=3' ]
FFfiles = [
'../2Mnacl_opc_Uext0.5_nacl_298K_Gauss/nacl_ff.dat',
'../2Mnacl_opc_Uext1_nacl_298K_Gauss/nacl_ff_1.dat',
'../2Mnacl_opc_Uext2_nacl_298K_Gauss/nacl_ff.dat',
'../2Mnacl_opc_Uext2.5_nacl_298K_Gauss/nacl_ff_1.dat',
'../2Mnacl_opc_Uext3_nacl_298K_Gauss/2Mnacl_opc_Uext3_nacl_298K_SplineToGauss_ff.dat']

PotName = 'LJGaussCl-_Na+'
PotType = 'gauss' #'gauss'
cut = 7.
rmax = cut

plotTitle = PotName
ylabel = 'Pair potential'
xlabel = '$r (\sigma)$'

##################
#check
if len(legends) != len(FFfiles):
    Exception('Mismatch in sizes of legends and FFfiles')
cwd = os.getcwd()

def getSpline(ParamDict, cut, rmax = 10, N = 1000, ParamName = 'Knots'):
    import spline
    us = []
    rs = np.linspace(0,rmax,N)
    if len(ParamDict.keys()) > 1:
        Exception('More than one potentials for spline')
    for Pot in ParamDict.keys():
        knots = ParamDict[Pot][ParamName]
        s = spline.Spline(cut, knots)
        for r in rs:
            us.append(s.Val(r))
    return rs, us
def getGauss(ParamDict, rmax = 10, N = 1000):
    rs = np.linspace(0,rmax,N)
    us = np.zeros(N)
    for Pot in ParamDict.keys():
        B = ParamDict[Pot]['B']
        Kappa = ParamDict[Pot]['Kappa']
        Dist0 = ParamDict[Pot]['Dist0']
        us += B * np.exp(-(rs-Dist0)**2 * Kappa)
    return rs, us

#get parameters from each ff file
ParamDicts = {}
xs = []
ys = []
for i, FFfile in enumerate(FFfiles): 
    num = []
    f = open(FFfile, 'r')
    str = f.read()
    f.close()
    ParamDict = {} #dictionary of paramters
    potential_str = [val for val in str.split('>>> POTENTIAL ') if len(val) > 0]
    for potential in potential_str:
        s = potential.split('{') #split potential name from params
        if  PotName in s[0]:
            PotNameTemp = s[0].split()[0]
            s = s[-1] #only take the parameters
            params = [s] #[' '.join(val.rsplit()) for val in re.split('}|,',s) if len(val)>0] 
            for i, param in enumerate(params):
                dict_val = ast.literal_eval('{'+param)
                ParamDict.update({PotNameTemp : dict_val})
    
    #each x value has a dictionary of parameters
    ParamDicts.update({legends[i]:ParamDict})

    if PotType == 'spline':
        rs, us = getSpline(ParamDict, cut, rmax = rmax)
    elif PotType == 'gauss':
        rs, us = getGauss(ParamDict, rmax = rmax)
    xs.append(rs)
    ys.append(us)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
for i, y in enumerate(ys):
    x = xs[i]
    ax.plot(x, y, marker='None', ls='-', lw=1, ms=4, label = legends[i])
plt.ylabel(ylabel)
plt.xlabel(xlabel)
title = plotTitle
plt.legend(loc='best', prop={'size': 5})
plt.title(title,loc='center')
plt.savefig('_'.join(title.split())+'.png',dpi=500,bbox_inches='tight')
plt.show()
    
