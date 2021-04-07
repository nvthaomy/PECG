import os, sys, re
import numpy as np
import mdtraj as md
import matplotlib
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
colors = ['#6495ED','r','#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#949c2d', '#a34a17','#c43b99','#949c2d','#1E90FF']

######

dirs = [
'3wtpercentPE_N30_f1_0.5MNaCl',
'3wtpercentPE_N30_f1_1MNaCl',
'3wtpercentPE_N30_f1_2MNaCl',
'3wtpercentPE_N30_f1_3MNaCl',
'3wtpercentPE_N30_f1_4MNaCl']
#'6wtpercentPE_N6_f1_0.5MNaCl',
#'6wtpercentPE_N12_f1_0.5MNaCl',
#'6wtpercentPE_N20_f1_0.5MNaCl',
#'6wtpercentPE_N30_f1_0.5MNaCl',
#'6wtpercentPE_N40_f1_0.5MNaCl']

#legends = ['N=6', 'N=12', 'N=20', 'N=30','N=40']
legends = ['0.5M','1M','2M','3M','4M']
pair = 'A- B+'
ext = 'N30 3wtpercentPE '
ylabel = 'RDF'

dataFiles = 'rdf_A-_B+.dat'
cwd = os.getcwd()
aw = 0.31 #nm
rs = []
gs = []

for dir in dirs:
    file = os.path.join(cwd,dir,dataFiles)
    rs.append(np.loadtxt(file)[40:,0]*aw)
    gs.append(np.loadtxt(file)[40:,1])

fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
axs.set_prop_cycle('color', colors)
#for i,r in enumerate(AArs):
#axs.plot(AArs[0], AAgs[0], marker=None, ls = '--', lw = 0.5, c= 'k', label = '4M AA')
#axs.plot(AArs[1], AAgs[1], marker=None, ls = '--', lw = 0.5, c= 'g', label = '4M AA')

for i,r in enumerate(rs):
    axs.plot(r, gs[i], marker=None,ls='-',lw=1.2, label = legends[i])
#plt.xlim((np.ravel(rs)).min(), 8.9) #(np.ravel(rs)).max())
plt.xlim(0.3,3)
plt.ylim(0)
plt.xlabel('r (nm)')
plt.ylabel(ylabel)
plt.legend(loc='best',prop={'size':5})
title = '{} rdf {}'.format(ext,pair)
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
