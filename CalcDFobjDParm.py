import os
import numpy as np
import scipy as sp
# conda install -c omnia mdtraj
import mdtraj as md


''' USER INPUTS '''
traj_cg = ['tmp04zVYP1.lammpstrj'] # the cg traj file    
  
traj_ref = 'testsim_trj.lammpstrj' # the reference traj file

cgTrajWarmup = 0
cgTrajSlice = 1
refTrajWarmup = 0
refTrajSlice = 1
ScaleRefTraj = 10. # incase you need to scale traj by 10. to get into Angstroms 
ScaleCGTraj  = 10. 
ShiftCoords = 0. # shift refTraj and CGTraj, not necessary if using .dcd files

# for custom topology
NPol = 1
DOP  = 20


# cg model bond parameters
FConst     = [2.0] # kbT/Ang**2
BondOffset = [2.5] # Ang
BondParam = 'Ro' # parameter to vary
RgTar = 8. # target Rg, Ang 
coef = 1. # for bias 


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

bondedpairs = np.asarray(bondedpairs,dtype='int32')
np.savetxt("bondedpairs.txt",bondedpairs)

def CalcDist(Pos1,Pos2,BoxL):
    '''
        Calculate the minimum distance under PBC 
    '''
    r1 = np.mod(Pos1,BoxL) 
    r2 = np.mod(Pos2,BoxL)
    r = r2 - r1

    r = r - BoxL*np.round(r/BoxL)
    return np.sqrt(np.sum(r**2))
    
def CalcDUBondDParam(_r,_k,_Ro,_nbonds,_param='Ro'):
    ''' 
        Returns the derivative of Ubond = k*(_r-_Ro)**2 with respect to _Ro or _k.
        
        _param = 'Ro' or 'k'
        
        DUbondDRo = -2*k*(_r-_Ro)
        
    '''
    _b = np.subtract(_r,_Ro)
    
    if _param == 'Ro':
        DUbondDRo = -2.*_k*_b
        DUbondDParam = DUbondDRo
    elif _param == 'k':
        _bsq = _b**2
        DUbondDk = _bsq
        DUbondDParam = DUbondDk
    else:
        print('BondParam {} not recognized.'.format(_param))
        
    return DUbondDParam
       
def CalculateRg(_traj,_scale):
    ''' 
        Calculate the Rg on each frame.
        Returns Avg Rg of each frame.
    '''

    RgData = []

    for _i,_chain in enumerate(_traj.topology.chains):
        _atoms = _traj.topology.select('chainid == {}'.format(_chain.index))
        Rg = md.compute_rg(_traj.atom_slice(_atoms))
        RgData.append(Rg*_scale)
        
    RgData = np.asarray(RgData)
    print('RgData Shape: {}'.format(RgData.shape))
    
    RgAvgPerFrame = np.average(RgData,axis=0)
    print('RgAvgPerFrame Shape: {}'.format(RgAvgPerFrame.shape))
    
    return RgAvgPerFrame

#Not Used Currently
def CalculateRgMinImage(_traj,_scale):
    ''' 
        Calculate the Rg on each frame.
        Returns Avg Rg of each frame.
    '''

    RgData = []
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
                # start from atom 1
                #atompos[:,_j,:] = _traj.xyz[:,_atoms[_j],:]
                continue
            else: 
                _p1 = _traj.xyz[:,int(_atomindx-1)]
                _p2 = _traj.xyz[:,int(_atomindx)]
                dPos = _p2 - _p1
                dPos = dPos - np.rint(dPos/_boxL)*_boxL
                atompos[:,int(_j),:] = atompos[:,int(_j-1),:] + dPos
        
        np.savetxt('atompos.txt',(atompos[:,0,0],atompos[:,1,0],atompos[:,2,0]))
        np.savetxt('atompos_traj.txt',_traj.xyz[0])
        
        _comperframe = np.average(atompos,axis=1) 
        print('COM Shape')
        print(_comperframe.shape)
        
        _dist = atompos - _comperframe[:,None,:]
        _distsq = np.sum(_dist*_dist,axis=2)
        
        print('Dist Sq Shape')
        print(_distsq.shape)
        
        _Rg = np.sqrt(np.sum(_distsq,axis=1)/len(_atoms))*_scale
        RgData[:,_i] = _Rg
        
        
    print('RgData Shape: {}'.format(RgData.shape))
    
    RgAvgPerFrame = np.average(RgData,axis=1)
    print('RgAvgPerFrame Shape: {}'.format(RgAvgPerFrame.shape))
    
    return RgAvgPerFrame

def CalcDRgDLambda(_DUDLambda,_RgAvgPerFrame):
    '''
        Calculate:
        Lambda = some CG FF parameter being optimized, e.g., Ro or FConst
        DRgDRo = <DU/DLambda>*<Rg> - <Rg*DU/DLambda>
    
    '''
    DUDLambdaRg = np.average((_DUDLambda*_RgAvgPerFrame))
    DRgDLambda = (np.average(_DUDLambda)*np.average(_RgAvgPerFrame)) - DUDLambdaRg

    return DRgDLambda,DUDLambdaRg
    
def CalcDFBiasDLambda(_coef,_RgAvg,_Rg_Tar,_DRgDLambdaAvg):
    '''
        Calculate DFbiasDLambda
        DFbiasDLambda = ceof*(<Rg>-Rg,T)*DRgDLambda
    
    '''
    
    DFbiasDLambda = _coef*(_RgAvg-_Rg_Tar)*_DRgDLambdaAvg
    
    return DFbiasDLambda
    
''' load in the ref trajectory '''
refTraj = md.load(traj_ref,top=top_,stride=int(refTrajSlice))
refTraj = refTraj[refTrajWarmup:]
refTraj.xyz[:] = refTraj.xyz[:]+ShiftCoords

#refTraj.save_lammpstrj('test.lammpstrj')
#np.savetxt('refTrajUnitcellvectors.txt',refTraj.unitcell_lengths)

print ("***Reference Trajectory***")
print ("Unit cell:")
print ("	{}".format(refTraj.unitcell_lengths[0])) # List of the unit cell on each frame
print ("Unit Cell Vectors:")
print ("	{}".format(refTraj.unitcell_vectors[0]))
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
DSrelDLambdaData = []
for _cgindex,_cgtraj in enumerate(traj_cg):
    ''' load in the trajectories '''
    cgTraj = md.load(_cgtraj,top=top_,stride=int(cgTrajSlice))
    cgTraj = cgTraj[cgTrajWarmup:]
    cgTraj.xyz[:] = cgTraj.xyz[:]+ShiftCoords
    
    
    ''' mdtraj analysis '''
    print ()
    print ("***Coarse-grain Trajectory # {}***".format(_cgindex))
    print ("Unit cell:")
    print ("	{}".format(cgTraj.unitcell_lengths[0])) # List of the unit cell on each frame
    print ("Unit Cell Vectors:")
    print ("	{}".format(cgTraj.unitcell_vectors[0]))
    print ('Number of frames:')
    print ("	{}".format(cgTraj.n_frames))
    print ('Number of molecule types:')
    print ("	{}".format(cgTraj.n_chains))
    print ('Number of molecules:')
    print ("	{}".format(cgTraj.n_residues))
    print ('Number of atoms:')
    print ("	{}".format(cgTraj.n_atoms))
    print ("Atom 1 coordinates:")
    print ('	{}'.format(cgTraj.xyz[0][0]))
        
        
    #Get distances for bonded pairs
    
    cgdist = md.compute_distances(cgTraj,bondedpairs,periodic=True,opt=True)
    np.savetxt('cgtraj_{}_dist_info.txt'.format(_cgindex),cgdist)  
    
    #Calculate the derivative of the Bond Potential in ref and cg model
    cgDUbond = np.sum(CalcDUBondDParam(cgdist*ScaleCGTraj,FConst[_cgindex],BondOffset[_cgindex],n_bonds,_param=BondParam),axis=1) # in units Angstrom
    refDUbond = np.sum(CalcDUBondDParam(refdist*ScaleRefTraj,FConst[_cgindex],BondOffset[_cgindex],n_bonds,_param=BondParam),axis=1) # in units Angstrom
    
    # Calculate <DUbondDLambda>_cgtraj
    cgDUbondAvg = np.average(cgDUbond)
    cgDUbondVar = np.var(cgDUbond)
    cgDUbondStdErr = np.sqrt(np.var(cgDUbond)/cgDUbond.shape[0])
    
    # Calculate <DUbondDLambda>_reftraj
    refDUbondAvg = np.average(refDUbond)
    refDUbondVar = np.var(cgDUbond)
    refDUbondStdErr = np.sqrt(np.var(refDUbond)/refDUbond.shape[0])
    
    # Calculate DSrelDLambda = <DUbondDLambda>_reftraj - <DUbondDLambda>_cgtraj
    DSrelDLambda = refDUbondAvg - cgDUbondAvg
    DSrelDRoStdErr = np.sqrt(cgDUbondStdErr**2+refDUbondStdErr**2)
    
    # Collect data vs Ro
    DSrelDLambdaData = [BondOffset[_cgindex],DSrelDLambda,DSrelDRoStdErr]

    DSrelDLambdaData = np.asarray(DSrelDLambdaData)
    np.savetxt('DSrelD{}.txt'.format(BondParam),DSrelDLambdaData)

    ''' Plot and Calculate DFobj/DRo & DFbias/DRo & DRgDRo '''

    # Calculate DRgDRo
    refRgAvgPerFrame = CalculateRg(refTraj,ScaleRefTraj)
    np.savetxt('refRgAvgPerFrame.txt',refRgAvgPerFrame)
    refRgAvg = np.average(refRgAvgPerFrame)
    
    RgAvgPerFrame = CalculateRg(cgTraj,ScaleCGTraj) 
    RgAvg = np.average(RgAvgPerFrame)
    
    DRgDLambdaAvg,DUDLambdaRg = CalcDRgDLambda(cgDUbond,RgAvgPerFrame)
    
    DFBiasDLambda = CalcDFBiasDLambda(coef,np.average(RgAvgPerFrame),RgTar,DRgDLambdaAvg)
    
    DFObjDLambda = DFBiasDLambda + DSrelDLambdaData[1]
    
    print('Calculations and Derivatives for Lambda {}:'.format(BondOffset[_cgindex]))
    print('     ref. RgAvg:         {}'.format(refRgAvg))
    print('     CG   RgAvg:         {}'.format(RgAvg))
    print('     DUD{}:         {}'.format(BondParam,cgDUbondAvg))
    print('     DUD{}*Rg:      {}'.format(BondParam,DUDLambdaRg))
    print('     DRgD{}:        {}'.format(BondParam,DRgDLambdaAvg))
    print('     DFBiasD{}:     {}'.format(BondParam,DFBiasDLambda))
    print('     DSrelD{}:      {}'.format(BondParam,DSrelDLambdaData[1]))
    print('     DFObjD{}:      {}'.format(BondParam,DFObjDLambda))
    
    with open('manualcalcs_values.txt','w') as fout:
        fout.write('    PBond.Dist0:    {}\n'.format(BondOffset))
        fout.write('    ref. RgAvg:     {}\n'.format(refRgAvg))
        fout.write('    CG.  RgAvg:     {}\n'.format(RgAvg))
        fout.write('    Rg  Target:     {}\n'.format(RgTar))
        fout.write('    Pen. Coef:      {}\n'.format(coef))
        fout.write('    DUD{}:          {}\n'.format(BondParam,cgDUbondAvg))
        fout.write('    DU{}*Rg:        {}\n'.format(BondParam,DUDLambdaRg))
        fout.write('    DRgD{}:         {}\n'.format(BondParam,DRgDLambdaAvg))
        fout.write('    DFBiasD{}:      {}\n'.format(BondParam,DFBiasDLambda))
        fout.write('    DSrelD{}:       {}\n'.format(BondParam,DSrelDLambdaData[1]))
        fout.write('    DUD{} refTraj:  {}\n'.format(BondParam,refDUbondAvg))
        fout.write('    DFObjD{}:       {}\n'.format(BondParam,DFObjDLambda))
    

