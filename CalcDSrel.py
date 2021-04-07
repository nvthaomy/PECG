import os
import numpy as np
import scipy as sp
import mdtraj as md
import re


''' USER INPUTS '''
varyFConst = False
varyRo = True
BondOffset_str = ['1.158', '1.20', '1.25', '1.30', '1.35', '1.4', '1.5', '1.6', '1.7', '1.8', '2.0', '2.2', '2.4', '2.6', '2.8', '3.0']
BondOffset = [float(f) for f in BondOffset_str]
traj_cg = []
for i in range(len(BondOffset)):
    traj_cg.append('../R0{}/PE0_traj.dcd'.format(BondOffset_str[i]))


traj_ref ='/home/mnguyen/PE/PAA_PAH/AA/xp0.1_10AA24f1_10AH24f1_325nacl_12500hoh/atacticV2/sim/TrajMap_traj_combined_200fr/TrajCOM.dcd'
# '/home/mnguyen/PE/PAA_PAH/CG/xp0.1_10AA24f1_10AH24f1_325nacl_12500hoh/1.200fCG_200fAA_cut10/atacticV2/2.400nsMD_fixedHOH_NaCl/trajectory_200nsMD_200frames.dcd'

rg_datafile = None
cg_der = None
cgTrajWarmup = 0
cgTrajSlice = 1
refTrajWarmup = 0
refTrajSlice = 5

NPol = 20
DOP  = 24
NWaters = 325*2 + 12500
Rgchainid = range(0,10)
_L =  235675.980541**(1./3.)
RgTar = 1.34283/0.31

# cg model bond parameters
FConst     = 5.1972e+01
BondOffset = BondOffset
# length scale: CG_length = scale * mdtraj_length
refscale = 1.0 #1/0.31
cgscale = 1.0

SimFobj_vs_Ro_Lag_Off = None #'SimFobj_vs_Ro_LagOff.txt'
SimFobj_vs_Ro_Lag_On = None #'SimFobj_vs_k_LagOn.txt'

#coefficients for penalty with Lagragian multiplier =  0
coef = [0,1,2,10,20,40,100,200,1000,1e4,1e5,1e6,1e7,1e8,1e9,1e10]
''' Build topology '''

top_ = md.Topology()
bondedpairs = []
cnt = 0
for i in range(NPol):
    top_.add_chain()
    top_.add_residue("POL",top_.chain(i))
    for j in range(DOP):
        top_.add_atom('P',md.element.carbon,top_.residue(i))
        if j > 0: # add bonds
            top_.add_bond(top_.atom(cnt-1),top_.atom(cnt))
            bondedpairs.append([top_.atom(cnt-1).index,top_.atom(cnt).index])
        cnt += 1

for i in range(NPol,NPol+NWaters):
    top_.add_chain()
    top_.add_residue("SOL",top_.chain(i))
    top_.add_atom('S',md.element.carbon,top_.residue(i))

bondedpairs = np.asarray(bondedpairs,dtype='int64')
np.savetxt("bondedpairs.txt",bondedpairs)

def CalcDRgDLambda(_DUDLambda,_RgAvgPerFrame):
    '''
        Calculate:
        Lambda = some CG FF parameter being optimized, e.g., Ro or FConst
        DRgDRo = <DU/DLambda>*<Rg> - <Rg*DU/DLambda>
    
    '''
    DUDLambdaRg = np.average((_DUDLambda*_RgAvgPerFrame))
    DRgDLambda = (np.average(_DUDLambda)*np.average(_RgAvgPerFrame)) - DUDLambdaRg

    return DRgDLambda,DUDLambdaRg
    
def CalcDFBiasDLambda(_coef,_lag,_RgAvg,_Rg_Tar,_DRgDLambdaAvg):
    '''
        Calculate DFbiasDLambda
        DFbiasDLambda = ceof*(<Rg>-Rg,T)*DRgDLambda
    
    '''
    
    DFbiasDLambda = _coef*(_RgAvg-_Rg_Tar)*_DRgDLambdaAvg - _lag * _DRgDLambdaAvg
    
    return DFbiasDLambda


def HarmonicBondDer(_r,_k,_Ro,_nbonds):
    ''' 
        Returns the derivative of Ubond = k*(_r-_Ro)**2 with respect to _Ro.
        
        DUbondDRo = -2*k*(_r-_Ro)
        
    '''
    _b = np.subtract(_r,_Ro)
    DUbondDRo = -2.*_k*_b    
    return DUbondDRo

def HarmonicBondDer2(_r,_k,_Ro,_nbonds):
    ''' 
        Returns the derivative of Ubond = k*(_r-_Ro)**2 with respect to _k.
        
        DUbondDRo = (_r-_Ro)^2
        
    '''
    _b = np.subtract(_r,_Ro)
    DUbondDk = _b**2    
    return DUbondDk

def CalculateRgMinImage(_traj,_scale,chainid=None):
    ''' 
        Calculate the Rg on each frame.
        Returns Avg Rg of each frame.
    '''

    if not chainid:
        _atoms = _traj.topology.select('chainid 0 to {}'.format(NPol-1))
    else:
        s = 'chainid '
        for _i in chainid:
            s += '{} '.format(_i)
        _atoms = _traj.topology.select(s)
    _traj = _traj.atom_slice(_atoms)
    RgData = np.zeros((_traj.n_frames,_traj.n_chains))
    
    _boxL  = _traj.unitcell_lengths
    np.savetxt('boxunitcell.txt',_boxL)
    for _i,_chain in enumerate(_traj.topology.chains):
        # get atoms in chain
        _atoms = _traj.topology.select('chainid == {}'.format(_chain.index))
        # initialize memory for min. imaged positions along chain backbone
        atompos = np.zeros((_traj.n_frames,len(_atoms),3))
        for _j,_atomindx in enumerate(_atoms):
            if _j == 0:
                continue
            else: 
                _p1 = _traj.xyz[:,int(_atomindx-1)]
                _p2 = _traj.xyz[:,int(_atomindx)]
                dPos = _p2 - _p1
                dPos = dPos - np.rint(dPos/_boxL)*_boxL
                atompos[:,int(_j),:] = atompos[:,int(_j-1),:] + dPos
                
                magdPos = np.sum(dPos*dPos,axis=1)
                magdPos = np.sqrt(magdPos)
        
        _comperframe = np.average(atompos,axis=1) 
        
        _dist = atompos - _comperframe[:,None,:]
        _distsq = np.sum(_dist*_dist,axis=2)
        
        _Rg = np.sqrt(np.sum(_distsq,axis=1)/len(_atoms))*_scale
        RgData[:,_i] = _Rg
    RgAvgPerFrame = np.average(RgData,axis=1)
    RgData = np.asarray(RgData)

    #print('RgData Shape: {}'.format(RgData.shape))
    #print('RgAvgPerFrame Shape: {}'.format(RgAvgPerFrame.shape))    
    return RgAvgPerFrame

''' load in the ref trajectory '''
refTraj = md.load(traj_ref,top=top_,stride=int(refTrajSlice))
refTraj = refTraj[refTrajWarmup:]
#RgAvgPerFrame = CalculateRgMinImage(refTraj,scale) 
RgAvgPerFrame = CalculateRgMinImage(refTraj,refscale,chainid=Rgchainid) 
RgAvg = np.average(RgAvgPerFrame) 
print('Reference Rg {}'.format(RgAvg))
 
print ("***Reference Trajectory***")
print ("Unit cell:")
print ("	{}".format(refTraj.unitcell_lengths[0])) # List of the unit cell on each frame
print ('Number of frames:')
print ("	{}".format(refTraj.n_frames))
print ('Number of molecule types:')
print ("	{}".format(refTraj.n_chains))
print ('Number of molecules:')
print ("	{}".format(refTraj.n_residues))
print ('Number of atoms:')
print ("	{}".format(refTraj.n_atoms))
print ("Atom 1 coordinates:")
print ('	{}'.format(refTraj.xyz[0][0]))
print ("Number of bonds")
n_bonds = len(bondedpairs)
print ('	{}'.format(n_bonds))

refdist = md.compute_distances(refTraj,bondedpairs,periodic=True,opt=True)
np.savetxt('ref_dist_info.txt',refdist)   

# iterate over cg_trajectories
DSrelDRoData = []
RgData = []
if cg_der:
    SimData = np.loadtxt(cg_der,delimiter=',') 
    coefSim = SimData[0,6]
    lagSim = SimData[0,7]
else:
    coefSim = 1000.
    lagSim = 0.
print('Parameters to plot derivatives: c = {} lag_mult = {}'.format(coefSim,lagSim))

for _cgindex,_cgtraj in enumerate(traj_cg):
    ''' load in the trajectories '''
    cgTraj = md.load(_cgtraj,top=top_,stride=int(cgTrajSlice))
    cgTraj = cgTraj[cgTrajWarmup:]

    ''' Put in unitcell info if there is none '''
    if cgTraj.unitcell_lengths is None:
        nframes = cgTraj.n_frames
        Lengths = np.asarray([[_L,_L,_L]]*nframes)
        cgTraj.unitcell_lengths = Lengths
        Volumes = np.asarray([[_L**3]]*nframes)
        #traj.unitcell_volumes = Volumes
    else:
        Lengths = cgTraj.unitcell_lengths
 
    if cgTraj.unitcell_vectors is None:
        nframes = cgTraj.n_frames
        Vectors = np.asarray([[[_L,0.,0.],[0.,_L,0.],[0.,0.,_L]]]*nframes)
        cgTraj.unitcell_vectors = Vectors
    #if traj.unitcell_angles is None:
    nframes = cgTraj.n_frames
    Angles = np.asarray([[90.,90.,90.]]*nframes)
        #traj.unitcell_angles = Angles
        
    traj = md.Trajectory(cgTraj.xyz, cgTraj.topology, time=None, unitcell_lengths=Lengths, unitcell_angles=Angles)

    cgTraj = traj

    ''' mdtraj analysis '''
    print ()
    print ("***Coarse-grain Trajectory # {}***".format(_cgindex))
    # print ("Unit cell:")
    # print ("	{}".format(cgTraj.unitcell_lengths[0])) # List of the unit cell on each frame
    # print ('Number of frames:')
    # print ("	{}".format(cgTraj.n_frames))
    # print ('Number of molecule types:')
    # print ("	{}".format(cgTraj.n_chains))
    # print ('Number of molecules:')
    # print ("	{}".format(cgTraj.n_residues))
    # print ('Number of atoms:')
    # print ("	{}".format(cgTraj.n_atoms))
    # print ("Atom 1 coordinates:")
    # print ('	{}'.format(cgTraj.xyz[0][0]))
        
        
    #Get distances for bonded pairs
    cgdist = md.compute_distances(cgTraj,bondedpairs,periodic=True,opt=True)
    np.savetxt('cgtraj_{}_dist_info.txt'.format(_cgindex),cgdist)  

    RgAvgPerFrame = CalculateRgMinImage(cgTraj,cgscale,chainid=Rgchainid) 
    RgAvg = np.average(RgAvgPerFrame) 
    print('Avg Rg: {}'.format(RgAvg))
    if varyRo:       
        RgData.append([BondOffset[_cgindex],RgAvg])
    elif varyFConst:
        RgData.append([FConst[_cgindex],RgAvg])
            
    if varyRo:
        #Calculate the derivative of the Bond Potential in ref and cg model
        cgDUbond = np.sum(HarmonicBondDer(cgdist*cgscale,FConst,BondOffset[_cgindex],n_bonds),axis=1) # in units Angstrom
        refDUbond = np.sum(HarmonicBondDer(refdist*refscale,FConst,BondOffset[_cgindex],n_bonds),axis=1) # in units Angstrom        

    elif varyFConst:
        #Calculate the derivative of the Bond Potential in ref and cg model
        cgDUbond = np.sum(HarmonicBondDer2(cgdist*cgscale,FConst[_cgindex],BondOffset,n_bonds),axis=1) # in units Angstrom
        refDUbond = np.sum(HarmonicBondDer2(refdist*refscale,FConst[_cgindex],BondOffset,n_bonds),axis=1) # in units Angstrom
        
    # Calculate <DUbondDRo>_cgtraj
    cgDUbondAvg = np.average(cgDUbond)
    cgDUbondVar = np.var(cgDUbond)
    cgDUbondStdErr = np.sqrt(np.var(cgDUbond)/cgDUbond.shape[0])
    
    # Calculate <DUbondDRo>_reftraj
    refDUbondAvg = np.average(refDUbond)
    refDUbondVar = np.var(cgDUbond)
    refDUbondStdErr = np.sqrt(np.var(refDUbond)/refDUbond.shape[0])
    
    # Calculate DSrelDRo = <DUbondDRo>_reftraj - <DUbondDRo>_cgtraj
    DSrelDRo = refDUbondAvg - cgDUbondAvg
    DSrelDRoStdErr = np.sqrt(cgDUbondStdErr**2+refDUbondStdErr**2)
    
    DRgDLambdaAvg,DUDLambdaRg = CalcDRgDLambda(cgDUbond,RgAvgPerFrame)
    DFBiasDLambda = CalcDFBiasDLambda(coefSim,lagSim, np.average(RgAvgPerFrame),RgTar,DRgDLambdaAvg)
    DFObjDLambda = DFBiasDLambda + DSrelDRo
    
    # Collect data vs Ro
    if varyRo:
        DSrelDRoData.append([BondOffset[_cgindex],DSrelDRo,DSrelDRoStdErr,DFBiasDLambda,DFObjDLambda,DRgDLambdaAvg])
    elif varyFConst:
        DSrelDRoData.append([FConst[_cgindex],DSrelDRo,DSrelDRoStdErr,DFBiasDLambda,DFObjDLambda,DRgDLambdaAvg])

DSrelDRoData = np.asarray(DSrelDRoData)
RgData = np.asarray(RgData)
np.savetxt('AvgRg.dat',RgData,header='Param Rg')
    
if varyRo:          
    np.savetxt('DSrelDRo.txt',DSrelDRoData)
elif varyFConst:
    np.savetxt('DSrelDFConst.txt',DSrelDRoData)
    
''' Plot and Integrate Data '''

import matplotlib as mpl
mpl.use('Agg')
#with open('mpl.data','w') as mplout:
#    mplout.write(str(mpl.rcParams))
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['axes.edgecolor'] = 'black'
mpl.rcParams['figure.dpi'] = '300'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
mpl.rcParams['axes.labelsize'] = '20'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.formatter.use_mathtext'] = True
from scipy.interpolate import interp1d
from scipy.integrate import simps
import matplotlib.pyplot as plt

DSrelDRo_Int = interp1d(DSrelDRoData[:,0],DSrelDRoData[:,1],kind='quadratic')
Ro_data = np.linspace(np.min(DSrelDRoData[:,0]),np.max(DSrelDRoData[:,0]),10000)

plt.semilogy(DSrelDRoData[:,0], DSrelDRoData[:,1], "ko",label='chk')
plt.semilogy(Ro_data, DSrelDRo_Int(Ro_data), "r--")
if cg_der:
    plt.semilogy(SimData[:,0],SimData[:,1],'gs',mfc='None',label='Sim')
#plt.title('Ref. Polymer-Solvent', fontsize=18)
plt.legend()
#plt.ylim(1e-3) #,1e5)
#plt.xlim(0,3.5)
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.ylabel('$dS_{rel}/dR_o$',fontsize=18)
    plt.savefig('DSrelDRo_semilogy.png')
elif varyFConst:
    plt.xlabel('$k$',fontsize=18)
    plt.ylabel('$dS_{rel}/dk$',fontsize=18)
    plt.savefig('DSrelDFConst_semilogy.png')

plt.close()

plt.plot(DSrelDRoData[:,0], DSrelDRoData[:,1], "ko",label='chk')
plt.plot(Ro_data, DSrelDRo_Int(Ro_data), "r--")
if cg_der:
    plt.plot(SimData[:,0],SimData[:,1],'gs',mfc='None',label='Sim')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.title('Ref. Polymer-Solvent', fontsize=18)
plt.legend()
xmin, xmax, ymin, ymax = plt.axis()
if ymin*ymax < 0.:
    plt.hlines(0.,xmin,xmax,colors = 'b')
#plt.ylim(1e-3,1e5)
#plt.xlim(0,3.5)
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.ylabel('$dS_{rel}/dR_o$',fontsize=18)
    plt.savefig('DSrelDRo.png')
elif varyFConst:
    plt.xlabel('$k$',fontsize=18)
    plt.ylabel('$dS_{rel}/dk$',fontsize=18)
    plt.savefig('DSrelDFConst.png')    
plt.close()
# integrate DSrelDRo to get dSrel

SrelData = []
SrelData.append([Ro_data[0],0.])
for _ind,_rup in enumerate(Ro_data.tolist()):
    if _ind < 10: 
        pass
    else:
        _int = simps(DSrelDRo_Int(Ro_data[0:_ind]),x=Ro_data[0:_ind])
        SrelData.append([Ro_data[_ind],_int])

# test delta Srel
if varyFConst:
    Ro_Start = min(FConst)
    Ro_End   = max(FConst)
elif varyRo:
    Ro_Start = min(BondOffset)
    Ro_End   = max(BondOffset)

RoGrid = np.linspace(Ro_Start,Ro_End,1000)
_intTest = simps(DSrelDRo_Int(RoGrid),x=RoGrid)
if varyRo:
    print('Srel change from Ro {} to Ro {}'.format(Ro_Start,Ro_End))
elif varyFConst:
    print('Srel change from FConst {} to FConst {}'.format(Ro_Start,Ro_End))
print(_intTest) # sim package gave: -9.4061 for dSrel


# find Srel Minimum and shift curve to zero it
SrelData = np.asarray(SrelData)
SrelMin = min(SrelData[:,1])
SrelData[:,1] =  SrelData[:,1] - SrelMin
   
plt.semilogy(SrelData[:,0], SrelData[:,1], "k-")
#plt.semilogy(3.2589,0.,'rx')
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.savefig('Srel_vs_Ro.png')
elif varyFConst:
    plt.xlabel('$k$',fontsize=18)
    plt.savefig('Srel_vs_FConst.png')
plt.ylabel(r'$\delta S_{rel}$',fontsize=18)
plt.close()

# plot Srel and Rg vs Ro


# interpolate RgData
Rg_Int = interp1d(RgData[:,0],RgData[:,1],kind='quadratic')
Ro_data = np.linspace(DSrelDRoData[0,0],DSrelDRoData[-1,0],100000)

fig,ax1 = plt.subplots()
if varyRo:
    ax1.semilogy(SrelData[:,0]/cgscale, SrelData[:,1], "k-",label='$S_{rel}$')
elif varyFConst:
    ax1.semilogy(SrelData[:,0]*cgscale**2, SrelData[:,1], "k-",label='$S_{rel}$')
ax1.semilogy(0.21214,0.,'rx',label='$S_{rel} opt. sim.$')
ax2 = ax1.twinx()
if varyRo:
    ax2.plot(RgData[:,0]/cgscale,RgData[:,1],"ro",label='$R_{g}$')
    ax2.plot(Ro_data/cgscale,Rg_Int(Ro_data),"r-")
    ax1.set_xlabel('$R_o$ $[nm]$',fontsize=18)
elif varyFConst:
    ax2.plot(RgData[:,0] * cgscale**2,RgData[:,1],"ro",label='$R_{g}$')
    ax2.plot(Ro_data * cgscale**2,Rg_Int(Ro_data),"r-") 
    

ax1.set_ylabel(r'$\delta S_{rel}$',fontsize=18)
#ax1.set_ylim(0.0001,10000)
plt.legend(loc=0)
if varyRo:
    plt.savefig('SrelRg_vs_Ro.png')
elif varyFConst:
    plt.savefig('SrelRg_vs_FConst.png')
plt.close()

# generate modified Srel Obj function with harmonic penalty
import copy
Srel_Int = interp1d(SrelData[:,0],SrelData[:,1],kind='quadratic')

_obj  = [] # list of interp1d objects
_objmin = [] # list of the Ro that minimizes the objective


# import Sim Srel Lag Off
if SimFobj_vs_Ro_Lag_Off:
    simData = np.loadtxt(SimFobj_vs_Ro_Lag_Off,delimiter=',')
else:
    simData = None
    
shiftmin = False
for _i,_coef in enumerate(coef):
    obj_temp = Srel_Int(Ro_data) + 0.5*_coef*(Rg_Int(Ro_data)*cgscale-RgTar)**2
    if shiftmin:
        obj_min = min(obj_temp)
        obj_temp = obj_temp - obj_min
    _obj.append(interp1d(Ro_data,obj_temp,kind='quadratic')) # deepcopy 2 enusre not passing by reference
    _objmin.append([_coef,min(_obj[-1](Ro_data)),Ro_data[np.argmin(_obj[-1](Ro_data))]])

# save Ro_min 
_objmin = np.asarray(_objmin)

if varyRo:
    np.savetxt('Fobj_vs_Ro.txt',_objmin,header='coef minObj minParam')
elif varyFConst:
    np.savetxt('Fobj_vs_FConst.txt',_objmin,header='coef minObj minParam')
    
    
# now construct plot of Fobj with varying coef    
color=plt.cm.coolwarm(np.linspace(0,1,len(coef)))
#plt.semilogy(SrelData[:,0], SrelData[:,1], "k-",label='$0$')
#plt.axvline(simData[0,2],color='k',linestyle='--',linewidth=0.75)
for _i,obj_Int in enumerate(_obj):
    if _i % 2 ==0:
        plt.semilogy(Ro_data, obj_Int(Ro_data), "k-", color=color[_i],label='{0:1.2e}'.format(coef[_i]))
    else:
        plt.semilogy(Ro_data, obj_Int(Ro_data), "k-", color=color[_i])
    try:
        plt.semilogy(simData[_i,2], obj_Int(simData[_i,2]), "o", color=color[_i])
#        plt.axvline(simData[_i+1,2],color=color[_i+1],linestyle='--',linewidth=0.75)
    except:
        pass
        
plt.ylabel(r'$F_{obj}$',fontsize=18)
#plt.ylim(1e-3,1e7)
plt.legend()
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.savefig('Srel_vs_Ro_Coef_semilogy.png')
elif varyFConst:
    plt.xlabel('$FConst$',fontsize=18)
    plt.savefig('Srel_vs_FConst_Coef_semilogy.png')
plt.close()

color=plt.cm.rainbow(np.linspace(0,1,len(coef)))
plt.plot(SrelData[:,0], SrelData[:,1], "k-",label='$0$')
for _i,obj_Int in enumerate(_obj):
    plt.plot(Ro_data, obj_Int(Ro_data), "k-", color=color[_i],label='{0:1.2e}'.format(coef[_i]))


plt.ylabel(r'$F_{obj}$',fontsize=18)
#plt.ylim(1e-1,1e8)
plt.legend()
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.savefig('Srel_vs_Ro_Coef.png')
elif varyFConst:
    plt.xlabel('$FConst$',fontsize=18)
    plt.savefig('Srel_vs_FConst_Coef.png')
plt.close()

#=== import Sim Srel Lag On ===
if SimFobj_vs_Ro_Lag_On:
    simDataLagOn = np.loadtxt(SimFobj_vs_Ro_Lag_On,delimiter=',')
    # Obj with Lag on
    coef2 = simDataLagOn[:,0]
    lag = simDataLagOn[:,1]

    _obj2  = [] # list of interp1d objects
    _objmin2 = [] # list of the Ro that minimizes the objective
    for _i,_coef in enumerate(coef2):
        obj_temp = Srel_Int(Ro_data) + 0.5*_coef*(Rg_Int(Ro_data)*cgscale-RgTar)**2 - lag[_i]*(Rg_Int(Ro_data)*cgscale-RgTar)
        if shiftmin:
            obj_min = min(obj_temp)
            obj_temp = obj_temp - obj_min
        _obj2.append(interp1d(Ro_data,obj_temp,kind='quadratic')) # deepcopy 2 enusre not passing by reference
        _objmin2.append([_coef,min(_obj2[-1](Ro_data)),Ro_data[np.argmin(_obj2[-1](Ro_data))]])
    
    _objmin2 = np.asarray(_objmin2)
    if varyRo:
        np.savetxt('Fobj_vs_Ro_LagOn.txt',_objmin2)
    elif varyFConst:
        np.savetxt('Fobj_vs_FConst_LagOn.txt',_objmin2)

    # now construct plot of Fobj with varying coef    
    color=plt.cm.coolwarm(np.linspace(0,1,len(coef2)))
    for _i,obj_Int in enumerate(_obj2):
        if _i % 5 ==0:
            plt.semilogy(Ro_data, obj_Int(Ro_data), "k-", color=color[_i],label='{0:1.2e}'.format(coef2[_i]))
        else:
            plt.semilogy(Ro_data, obj_Int(Ro_data), "k-", color=color[_i])
        try:
            #plt.axvline(simDataLagOn[_i+1,2],color=color[_i+1],linestyle='--',linewidth=0.75)
            plt.semilogy(simDataLagOn[_i,2], obj_Int(simDataLagOn[_i,2]), "o", color=color[_i])
        except:
            pass
        
    plt.ylabel(r'$F_{obj}$',fontsize=18)
    #plt.ylim(1e-3,1e7)
    plt.legend(loc='lower left')
    if varyRo:
        plt.xlabel('$R_o$',fontsize=18)
        plt.savefig('Srel_vs_Ro_Coef_LagOn_semilogy.png')
    elif varyFConst:
        plt.xlabel('$FConst$',fontsize=18)
        plt.savefig('Srel_vs_FConst_Coef_LagOn_semilogy.png')
    plt.close()

    color=plt.cm.rainbow(np.linspace(0,1,len(coef2)))
    plt.plot(SrelData[:,0], SrelData[:,1], "k-",label='$0$')
    for _i,obj_Int in enumerate(_obj2):
        plt.plot(Ro_data, obj_Int(Ro_data), "k-", color=color[_i],label='{0:1.2e}'.format(coef2[_i]))
    plt.ylabel(r'$F_{obj}$',fontsize=18)
    #plt.ylim(1e-1,1e8)
    plt.legend()
    if varyRo:
        plt.xlabel('$R_o$',fontsize=18)
        plt.savefig('Srel_vs_Ro_Coef_LagOn.png')
    elif varyFConst:
        plt.xlabel('$FConst$',fontsize=18)
        plt.savefig('Srel_vs_FConst_Coef_LagOn.png')
    plt.close()

    
# plot Ro vs Coef 
try:
    plt.semilogx(simData[:,0],simData[:,2],"ko-",label='sim. Lag. = 0.')
except:
    pass
if SimFobj_vs_Ro_Lag_On:
    plt.semilogx(simDataLagOn[:,0],simDataLagOn[:,2],"g*-",label='sim. Lag. On')
    plt.semilogx(_objmin2[:,0],_objmin2[:,2],"bo-",label='chk Lag. On')
plt.semilogx(_objmin[:,0],_objmin[:,2],"rx-",label='chk Lag. Off')
plt.legend()
plt.xlabel('$C$')
if varyRo:
    plt.ylabel('$R_o$')
    plt.savefig('Ro_vs_Coef.png')
elif varyFConst:
    plt.ylabel('$FConst$')
    plt.savefig('FConst_vs_Coef.png')    
plt.close()

''' Calculate Derivative of Fobj '''
from scipy.misc import derivative as der
from scipy.interpolate import splrep as splineInt
from scipy.interpolate import splev

Romin = Ro_data[0]
Romax = Ro_data[-1]
dx = 1e-4
num = int((Romax-2*dx)/dx)
#Ro_new = np.linspace(Romin+dx,Romax-dx,num)
#x = Ro_new[1]-Ro_new[0]
print('dx: {}'.format(dx))

plt.plot(DSrelDRoData[:,0],DSrelDRoData[:,4],'ko',label='chk. Coef={0:1.1e} Lag.={1:1.1e}'.format(coefSim,lagSim))
if cg_der:
    plt.plot(SimData[:,0],SimData[:,2],'gs',mfc='None',label='Sim')
xmin, xmax, ymin, ymax = plt.axis()
if ymin*ymax < 0.:
    plt.hlines(0.,xmin,xmax,colors = 'b')

plt.legend()
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.ylabel(r'$dF_{obj}/dR_{o}$',fontsize=18)
    plt.savefig('DFobjDRo.png')
elif varyFConst:
    plt.xlabel('$FConst$',fontsize=18)
    plt.ylabel(r'$dF_{obj}/dFConst$',fontsize=18)
    plt.savefig('DFobjFConst.png')
plt.close()

plt.plot(DSrelDRoData[:,0],DSrelDRoData[:,3],'ko',label='chk. Coef={0:1.1e} Lag.={1:1.1e}'.format(coefSim,lagSim))
if cg_der:
    plt.plot(SimData[:,0],SimData[:,3],'gs',mfc='None',label='Sim')
xmin, xmax, ymin, ymax = plt.axis()
if ymin*ymax < 0.:
    plt.hlines(0.,xmin,xmax,colors = 'b')

plt.legend()
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.ylabel(r'$dF_{bias}/dR_{o}$',fontsize=18)
    plt.savefig('DFbiasDRo.png')
elif varyFConst:
    plt.xlabel('$FConst$',fontsize=18)
    plt.ylabel(r'$dF_{bias}/dFConst$',fontsize=18)
    plt.savefig('DFbiasFConst.png')
plt.close()

plt.plot(DSrelDRoData[:,0],DSrelDRoData[:,5],'ko',label='chk')
if cg_der:
    plt.plot(SimData[:,0],SimData[:,4],'gs',mfc='None',label='Sim')
xmin, xmax, ymin, ymax = plt.axis()
if ymin*ymax < 0.:
    plt.hlines(0.,xmin,xmax,colors = 'b')
plt.legend()
if varyRo:
    plt.xlabel('$R_o$',fontsize=18)
    plt.ylabel(r'$dRg/dR_{o}$',fontsize=18)
    plt.savefig('DRgDRo.png')
elif varyFConst:
    plt.xlabel('$FConst$',fontsize=18)
    plt.ylabel(r'$dRg/dFConst$',fontsize=18)
    plt.savefig('DRgFConst.png')





