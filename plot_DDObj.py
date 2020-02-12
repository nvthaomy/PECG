""" read DDObj from output file and plot its trajectory"""
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

plotTitle = 'DDObs vs Uext'
ylabel = 'DDObj'
xlabel = '$U_{ext}$'
xs = [0,0.5,1,2,2.5,3]
Dirs = [
'2Mnacl_opc_298K_NPT_fixUwater/2.ElecSys_200f_cut6/podsubmit.sh.o1103076',
'2Mnacl_opc_Uext0.5_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6/podsubmit.sh.o1105446',
'2Mnacl_opc_Uext1_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6/podsubmit.sh.o1103179',
'2Mnacl_opc_Uext2_NaCl_298K_NVT_fixUwater/7.ElecSys_200f_cut6/podsubmit.sh.o1103172',
'2Mnacl_opc_Uext2.5_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6/podsubmit.sh.o1103178',
'2Mnacl_opc_Uext3_NaCl_298K_NVT_fixUwater/2.ElecSys_200f_cut6/podsubmit.sh.o1103186']

dataFName = ''
keywords = 'DDObj'
#indices of iteration to plot
iterIds = [0,1,2,3,4]
##################
#check
if len(xs) != len(Dirs):
    Exception('Mismatch in sizes of x and Dirs')

reading =  False

#DDObjDicts = {Dir1_xlabel: {param1: [DDObj_Iter1, DDObj_Iter2, ...], param2: []}, Dir2_xlabel: {}} 
DDObjDicts = {}

cwd = os.getcwd()
for i, dir in enumerate(Dirs): 
    print('reading {}'.format(dir))
    num = []
#    DDObjDict = {param1: [DDObj_Iter1, DDObj_Iter2, ...], param2: []}
    DDObjDict = {}
#    f = open(os.path.join(dir+'/',dataFName), 'r')
    f = open(dir,'r')
    lines = f.readlines()
    for line in lines:
        vals = line.split()
        if keywords in line:
            colId = vals.index(keywords)
            reading = True
        if reading and len(vals) >= colId and not keywords in line: 
            param = vals[0]
            DDObj = float(vals[colId])
            if not param in DDObjDict.keys():
                DDObjDict.update({param: [DDObj]})
            else:
                DDObjDict[param].append(DDObj)
        elif reading and len(vals) == 0: #stop appending DDObj when reach an empty line
            reading = False

    DDObjDicts.update({xs[i]: DDObjDict})

print(DDObjDicts)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
for iterId in iterIds:
    ys = []
    xs_plot = []
    for i,x in enumerate(xs):
        DDObjDict = DDObjDicts[x]
        y = 0
        try:
             for param, DDObjs in DDObjDict.items(): 
                 y += DDObjs[iterId]
             xs_plot.append(x)
             ys.append(y)
        except:
             pass
    ax.plot(xs_plot, ys, marker='o', ls=':', lw=1, ms=4, label='IterId={}'.format(iterId))
plt.ylabel(ylabel)
plt.xlabel(xlabel)
plt.legend(loc='best')
title = plotTitle
plt.title(title,loc='center')
plt.savefig('_'.join(title.split())+'.png',dpi=500,bbox_inches='tight')
plt.show()
    
