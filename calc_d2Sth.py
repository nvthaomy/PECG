# total = thermo + residual part
# thermo = mean + susceptibility part = total - residual part
import numpy as np
import sys

label = sys.argv[1]
hessian_file = str(sys.argv[2]) #'nacl_recalc_Hessian_masked.dat'
hessian_relaxed_file = str(sys.argv[3]) #'relaxUext/nacl_recalc_Hessian_masked.dat'

d0 = np.loadtxt(hessian_file)
dr = np.loadtxt(hessian_relaxed_file)
#print('--- initial data ---\n{}'.format(d))

paraminds = np.array([0,1,2,3,4])
nparam = len(paraminds)

#in this case, only have one 'thermo parameter', we need to pair-up/collapse the last row and column into the 2nd-to-last row and column
#dnew = d[:-1,:]
#new[-1,:] += d[-1,:]
#dnew2 = dnew[:,:-1]
#dnew2[:,-1] += dnew[:,-1]
#d = dnew2
#print('--- processed data ---\n{}'.format(d))

indices = np.arange(d0.shape[0])
mask = np.ones(indices.shape, bool)
mask[paraminds]=False

param_indices = indices[np.invert(mask)]
thermo_indices = indices[mask]

def thermo_reduce( d2Srel0, d2Srel1 ):
    """ '0' is the actual state, '1' is the relaxed-Uext state """
    #thermo_matrix = d[thermo_indices[:,None],thermo_indices]
    thermo_matrix = d2Srel1[np.ix_(thermo_indices,thermo_indices)]

    d2Sth_mean = d2Srel0[np.ix_(param_indices,param_indices)] - d2Srel1[np.ix_(param_indices,param_indices)]
    d2Sth_fluct = np.zeros( ( nparam,nparam ) )
    for ii,ind1 in enumerate(paraminds):
        for jj,ind2 in enumerate(paraminds):
            v1 = d2Srel1[ind1,thermo_indices]
            v2 = d2Srel1[thermo_indices,ind2]
            tmp = np.matmul(np.linalg.inv(thermo_matrix),v2)
            d2Sth_fluct[ii,jj] = np.matmul(v1,tmp)

    print('--- d2Sth from macrostate error ---')
    print(d2Sth_mean)
    print('--- d2Sth from thermo susceptibility ---')
    print(d2Sth_fluct)

    d2Sth = d2Sth_mean + d2Sth_fluct
    return d2Sth, d2Sth_mean, d2Sth_fluct

d2Sth, d2Sth_m, d2Sth_chi = thermo_reduce(d0,dr)

if d2Sth.size == 1:
    print('total\tthermo')
    print('{}\t{}'.format(d0[0,0],d2Sth[0,0]))
else:
    print('--- eigenvalue analysis ---')
    d2Sres = d0[np.ix_(paraminds,paraminds)] - d2Sth
    et,dt = np.linalg.eig(d0[np.ix_(paraminds,paraminds)]) #total matrix eig
    er,dr = np.linalg.eig(d2Sres) #residual matrix eig
    eth,dth = np.linalg.eig(d2Sth) #thermo matrix eig

    indmax = np.argmax(eth)
    print('...eigenvalues\n{}'.format(eth))
    print('...eigenvectors\n{}'.format(dth))
    print('...max eigenvector: {}'.format(dth[:,indmax]))

    ratio = eth[indmax] / np.linalg.norm( d2Sres @ dth[:,indmax] )
    #print('{}\t{}'.format(ratio,np.max(et))
    print('max eigen value {}'.format(np.max(eth)))

    a =  np.linalg.norm( d2Sth @ dth[:,indmax] )
    b = np.linalg.norm( d2Sres @ dth[:,indmax] )
    print('np.linalg.norm( d2Sres @ dth[:,indmax] ): {}'.format(b))    
    print('thermo info/ residual info : {}'.format(ratio))
    print('thermo info/ residual info : {}'.format(a/b))

    np.savetxt('{}_d2S_tot.txt'.format(label),d0[np.ix_(paraminds,paraminds)])
    np.savetxt('{}_d2S_th.txt'.format(label),d2Sth)
    np.savetxt('{}_d2S_chi.txt'.format(label),d2Sth_chi)
    np.savetxt('{}_d2S_m.txt'.format(label),d2Sth_m)
    np.savetxt('{}_d2S_res.txt'.format(label),d2Sres)

    np.savetxt('{}_d2Sth_eigs.txt'.format(label),np.vstack([eth,dth]))

