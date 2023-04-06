''' compute RDF of models in the expanded ensemble and add them up using provided weights
'''
import os
from subprocess import call
import numpy as np 
import time
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)

prefix = 'target'
trajs = ['trajectory_100-825ns_3626fr_mapped.lammpstrj',
'trajectory_100-822ns_3609fr_mapped.lammpstrj',
'trajectory_100-733ns_3569fr_mapped.lammpstrj']

#prefix = 'model'
#trajs = ['../tmpmuRRKv1.trj.dcd',
#'../tmpSisvli1.trj.dcd',
#'../tmp3VbsVm1.trj.dcd']

tops = ['../tmp4U_pol1.initial.pdb',
'../tmpUwchy41.initial.pdb',
'../tmpCOHr_E1.initial.pdb']

weights = [1.,1.,1.]
ff_bead_map = {'A1': ['W','Y','I','L','V'], 'A2': ['A','G','P'], 'A3': ['S','E','K','Q'], 'HOH': ['HOH']}
atoms = [a for a in ff_bead_map.keys()]
print('Atom types for rdf calculation: ',atoms)
hist = {}
for k,traj in enumerate(trajs):
    top = tops[k]
    if not os.path.isfile(traj):
        call('gunzip {}'.format(traj + '.gz'), shell=True)
    for i in range(len(atoms)-1): #skip HOH-HOH
        a1 = ' '.join(ff_bead_map[atoms[i]])
        a1_ = '_'.join(ff_bead_map[atoms[i]])
        for j in range(i,len(atoms)): 
            a2 = ' '.join(ff_bead_map[atoms[j]])
            a2_ = '_'.join(ff_bead_map[atoms[j]])

            # remove file if already exists:
            if os.path.isfile('rdf_{}-{}.dat'.format(a1_,a2_)):
                os.remove('rdf_{}-{}.dat'.format(a1_,a2_))

            if 'HOH' in [atoms[i], atoms[j]]:
                s = 'python /home/mnguyen/bin/scripts/mdtraj_rdf.py {} {} -s 10 -a1 {} -a2 {} -c -rmax 10'.format(traj, top, a1, a2)
            else:
                s = 'python /home/mnguyen/bin/scripts/mdtraj_rdf.py {} {} -a1 {} -a2 {} -c -rmax 10'.format(traj, top, a1, a2)
            print(s)
            call(s, shell=True)
            while not os.path.exists('rdf_{}-{}.dat'.format(a1_,a2_)):
                time.sleep(1)
            os.rename('rdf_{}-{}.dat'.format(a1_,a2_), prefix + '{}_rdf_{}-{}.dat'.format(k,atoms[i],atoms[j]))
            os.rename('rdf_{}-{}.png'.format(a1_,a2_), prefix + '{}_rdf_{}-{}.png'.format(k,atoms[i],atoms[j]))
            os.rename('coord_{}-{}.png'.format(a1_,a2_), prefix + '{}_coord_{}-{}.png'.format(k,atoms[i],atoms[j]))
            rdf_fn = prefix + '{}_rdf_{}-{}.dat'.format(k,atoms[i],atoms[j])

            if k == 0:
                hist.update({(atoms[i],atoms[j]): [np.loadtxt(rdf_fn)[:,0], weights[i]/np.sum(weights) * np.loadtxt(rdf_fn)[:,1], weights[i]/np.sum(weights) * np.loadtxt(rdf_fn)[:,2]]}) # x, rdf, coord
            else: #update current histogram
                hist[(atoms[i],atoms[j])][1] += weights[i]/np.sum(weights) * np.loadtxt(rdf_fn)[:,1]
                hist[(atoms[i],atoms[j])][2] += weights[i]/np.sum(weights) * np.loadtxt(rdf_fn)[:,2]
            
            fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2]) 
            plt.plot(hist[(atoms[i],atoms[j])][0], hist[(atoms[i],atoms[j])][1], lw=1, c='k', label='total')    
            plt.xlabel('r')
            plt.ylabel('g(r)')
            plt.savefig(prefix+'_rdf_{}-{}.png'.format(atoms[i],atoms[j]),dpi=500,transparent=False,bbox_inches="tight")

            fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2]) 
            plt.plot(hist[(atoms[i],atoms[j])][0], hist[(atoms[i],atoms[j])][2], lw=1, c='k', label='total')    
            plt.xlabel('r')
            plt.ylabel('n(r)')
            plt.savefig(prefix+'_coord_{}-{}.png'.format(atoms[i],atoms[j]),dpi=500,transparent=False,bbox_inches="tight")

            if k == len(trajs) - 1:
                np.savetxt(prefix+'_total_rdf_{}-{}.dat'.format(atoms[i],atoms[j]), 
                            np.array(hist[(atoms[i],atoms[j])]).T, header='x rdf coord.')
