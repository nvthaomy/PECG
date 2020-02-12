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

plotTitle = 'LJGaussCl-_Na+0:B vs Uext'
ylabel = 'LJGaussCl-_HOH0:B'
xlabel = 'Iteration'
legends = [0.5,1,2,2.5,3]
Dirs = [
#'2Mnacl_opc_298K_NPT_fixUwater/2.ElecSys_200f_cut6',
'2Mnacl_opc_Uext0.5_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6',
'2Mnacl_opc_Uext1_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6',
'2Mnacl_opc_Uext2_NaCl_298K_NVT_fixUwater/7.ElecSys_200f_cut6',
'2Mnacl_opc_Uext2.5_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6',
'2Mnacl_opc_Uext3_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6']

dataFName = 'nacl_log.txt'
keywords = 'LJGaussCl-_Cl-0:B'
##################
#check
if len(legends) != len(Dirs):
    Exception('Mismatch in sizes of x and Dirs')

#ValDicts = {Dir1_xlabel: [Val_Iter1, Val_Iter2, ...],  Dir2_xlabel: [...]} 
ValDicts = {}

cwd = os.getcwd()
for i, dir in enumerate(Dirs): 
    Val = []
    f = open(os.path.join(dir+'/',dataFName), 'r')
#    f = open(dir,'r')
    lines = f.readlines()
    for j, line in enumerate(lines):
        vals = line.split()
        if j == 0:
            colId = vals.index(keywords)
        else: 
            Val.append(float(vals[colId]))

    ValDicts.update({legends[i]:Val})


fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
for i,legend in enumerate(legends):
    ys = ValDicts[legend]
    xs = range(len(ys))
    ax.plot(xs, ys, marker='None', ls='-', lw=1, ms=4, label=legend)
plt.ylabel(ylabel)
plt.xlabel(xlabel)
plt.legend(loc='best')
title = plotTitle
plt.title(title,loc='center')
plt.savefig('_'.join(title.split())+'.png',dpi=500,bbox_inches='tight')
plt.show()
    
