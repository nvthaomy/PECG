'''
plot histograms of targets and models outputted by sim
(average with equal weights if detect multiple targets in the case of extended ensemble)
'''

import numpy as np 
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)

mod_hist_file = open('poly_modhist.txt','r')
lines = mod_hist_file.readlines()
mod_hist_file.close()

hist = {}
read = False
for j,line in enumerate(lines):
    if line.startswith('POTENTIAL'):
        potential = line.split()[1]
        print(line)
    elif line.startswith('arg'):
        read = True
        x = []
        y = []
        z = []
    elif read:
        vals = [float(a) for a in line.split()]
        if len(line.split()) > 0:
            x.append(vals[0])
            y.append(vals[1::2]) # collect all targets
            z.append(vals[2::2]) # collect all models
            if len(hist) == 0 and len(x) == 1: print('detect {} targets'.format(len(vals[1::2])))
        if len(line.split()) == 0 or j == len(lines)-1:
            y = np.array(y) 
            z = np.array(z)
            # plot targets
            fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
            if y.shape[1] > 1:
                for i in range(y.shape[1]):
                    plt.plot(x,y[:,i]/y.shape[1], lw=0.75, ls='--')
                plt.plot(0,0, lw=0.75, ls='--', c='k', label='model * weight')
            y = np.mean(y, axis=1) #take average of all targets, assume equal weights 
            plt.plot(x,y, lw=1, c='k', label='total')    
            plt.legend(loc='upper right',prop={'size':5})
            plt.xlabel('x')
            plt.ylabel('distribution')
            plt.savefig(potential+'_tar.png',dpi=500,transparent=False,bbox_inches="tight")
            # plot models
            fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
            if z.shape[1] > 1:
                for i in range(z.shape[1]):
                    plt.plot(x,z[:,i]/z.shape[1], lw=0.75, ls='--')
                plt.plot(0,0, lw=0.75, ls='--', c='k', label='model * weight')
            z = np.mean(z, axis=1) #take average of all targets, assume equal weights 
            plt.plot(x,z, lw=1, c='k', label='total')    
            plt.legend(loc='upper right',prop={'size':5})
            plt.xlabel('x')
            plt.ylabel('distribution')
            plt.savefig(potential+'_mod.png',dpi=500,transparent=False,bbox_inches="tight")

            hist.update({potential : [x,y,z]})
            np.savetxt(potential+'.txt', np.array([x,y,z]).T, header = 'x target_hist model_hist')
            read = False
#plt.show()