#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:11:03 2019

@author: my
"""
import sys
sys.path.append('/home/mnguyen/bin/PECG/')
import sim, pickleTraj
import mappoly
import system
import optimizer
import numpy as np
"""
run with mappoly.py (or mappoly_implicit.py), system.py, optimizer.py, forcefield.py
Run expanded ensemble Srel opt. (number of System can be 1) with/without weights in NPT or NVT
If CGtrajs is empty, map AA traj to CG traj and scale positions and box size by 1/lengthScale (needed when running srel in dimensionless units)

Molecules:
    Na+, Cl-
    Polycation: B, B+
    Polyanion: A, A-
    water: HOH    
Supported potentials:
    Bond: harmonic bond with or without offset
    Pair: LJGauss (any number of gaussians) or spline
    Electrostatics: smeared Coulomb
    Sinusoidal external potential

****assumptions:
currently only do 1AAmonomer:1CGbead mapping for polymer
all atom indices in trajectory follow the same order as definition of MolTypes, NMons, and NMols
always fix gaussian strength for HOH-HOH
all pair types have the same fix
cross-a_ev = mean of self a_ev's
Born radius = a_el = sqrt(pi) * a_ev
if use multiple gaussians:
    Gaussians with even idices are repulsive, odd indices are attractive (by setting B.Min, B.Max)

****To do:
   allow mapping other than 1:1 for polymer
   read bond info with mdtraj to calculate number of monomers (NMons) ber bead and NMols
   how to run neutral system along with charged system? lammps doesn't allow coulomb for neutral beads
"""

"""---Inputs---"""  
kB = 1.380658e-23
TempRef = 298.15
kT = kB*TempRef
kTkJmol = kT/1000.0*6.022e23
kTkcalmol = kTkJmol/4.184
Units = sim.units.DimensionlessUnits #sim.units.AtomicUnits
lengthScale  = 3.1 # = CG length scale/ AA length scale = 3.1 A/ 1 A
if Units == sim.units.DimensionlessUnits:
    lengthScale = lengthScale
else:
    lengthScale = 1.
#dimensionless scales:
# length: a_water = 3.1 A
# energy: kT
# pressure: kT/a_water**3

"""TOPOLOGY"""
UseLJGauss_MolType = False
UseLJGauss_CustomBead = True
# define bead type for LJGauss
if UseLJGauss_CustomBead:
    ff_bead_map = {'A1': ['W','Y','I','L','V'], 'A2': ['A','G','P'], 'A3': ['S','E','K','Q'], 'HOH': ['HOH']}

#map AA traj to CG traj
AAtrajs = ['/home/mnguyen/BioSCOUT/openMM/4IDP-NME_2chains_0MNaCl_9nm/sim/analysis/trajectory_100-365ns_2650fr.dcd','/home/mnguyen/BioSCOUT/openMM/ACE-YWSA-2Yx2A_2chains_0MNaCl_9nm/sim/analysis/trajectory_100-521ns_4210fr.dcd','/home/mnguyen/BioSCOUT/openMM/4IDP-NME_ACE-YWSA-2Yx2A_9nm/sim/analysis/trajectory_100-317ns_2165fr.dcd']
AAtops = [ '/home/mnguyen/BioSCOUT/openMM/4IDP-NME_2chains_0MNaCl_9nm/sim/run1/top.pdb','/home/mnguyen/BioSCOUT/openMM/ACE-YWSA-2Yx2A_2chains_0MNaCl_9nm/sim/run1/top.pdb', '/home/mnguyen/BioSCOUT/openMM/4IDP-NME_ACE-YWSA-2Yx2A_9nm/sim/run1/top.pdb']
strides = [3, 4, 2]

#map AA residue to CG bead name
nameMap = {'Na+':'Na+', 'Cl-':'Cl-', 'HOH': 'HOH', 'WAT': 'HOH', 'NA+':'Na+', 'CL-':'Cl-',
           'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N', 'GLN': 'Q', 'LYS': 'K', 'THR': 'T', 'PRO': 'P', 'ASX': 'B', 'HIS': 'H', 'PHE': 'F', 'ALA': 'A', 'GLY': 'G', 'ILE': 'I', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'ACE': 'skip', 'NME': 'skip'}

#CG trajectories, if CGtrajs = [] will map AA traj, else will use provided CGtrajs
CGtrajs = []
#provide UniqueCGatomTypes if CGtrajs is not empty
UniqueCGatomTypes = []
# topology of molecules in all systems
head = ['S', 'P', 'A', 'E', 'A', 'K', 'S', 'P', 'V', 'E', 'V', 'K', 'S', 'P', 'A', 'E', 'A', 'K', 'S', 'P', 'V', 'E', 'V', 'K']
tail = ['Y', 'W', 'S', 'A', 'Y', 'G', 'A', 'Y', 'A', 'Q', 'Y', 'V', 'Y', 'I', 'Y', 'A', 'Y', 'W', 'Y', 'L']
UniqueMolTypes = {'Head': head, 'Tail': tail,'HOH':['HOH']}

#name of molecules in systems
#must in the right sequence as molecules in the trajectory
MolNamesList = [['Head','HOH'], ['Tail','HOH'], ['Head','Tail','HOH']]

# number of molecules for each molecule type, nSys x molecule types
NMolsDicts = [{'Head':2,'HOH': 25500}, {'Tail': 2, 'HOH': 25500}, {'Head': 1, 'Tail': 1, 'HOH': 25500}]

charges1 = {'Na+': 1., 'Cl-': -1., 'HOH': 0.}
amino_acids = ['A', 'C', 'B', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', 'X', 'Z', 'Z+', 'Z-']
charges2 = {key: 0. for key in amino_acids}
for key in ['R','H','K', 'Z+']:
    charges2[key] = 1.0
for key in ['D','E', 'Z-']:
    charges2[key] = -1.0
charges = {key: val for d in (charges1, charges2) for key, val in d.items()}
chargesNeutral  = {key: 0. for key in charges.keys()}

elements1 = {key: key for key in ['Na+','Cl-','HOH']}
elements2 = {key: key for key in amino_acids}
#elements2 = {key: 'Head' for key in amino_acids}
#elements2['X'] = 'Tail'
elements = {key: val for d in (elements1,elements2) for key, val in d.items()}
Name = 'poly'

"""INTEGRATION PARAMS"""
#real units: dt (ps), temp (K), pressure (atm)
dt = 0.05
#TempSet and PresSet are lists 
TempSet = [1.,1.,1.]
PresSet = [3.218,3.218,3.218] #enter values to enable NPT

PresSet = np.array(PresSet)
#PresSet *= 101325. #joules/m**3
#PresSet *= (10.**-10.)**3. * 0.000239006 *6.022e23 #kcal/mol/A**3

IntParams = {'TimeStep': dt, 'LangevinGamma': 1/(100*dt)}

"""SREL OPT"""
ConstrainMode = 1
SrelOnNeutralSys = True
UseLammps = False
UseOMM = True
UseSim = False
ScaleRuns = True
StepScales = [] #set to empty if don't want to scale steps
StepsEquil = int(1000/dt)
StepsProd = int(20000/dt)
StepsStride = int(StepsProd/500)
WeightSysByNMol = False
WeightSysByNAtom = False
RandomMul = 100
UseWPenalty = False
StageCoefs = [1e3,1e5, 1e8, 1e10]

MaxIter=None
SteepestIter=0

RgConstrain=False
MolIdRgList=[[range(0,10),range(10,20)],[],[]]
RgTarList = [[1.3568/0.31,1.1892/0.31],[],[]]

"""FORCEFIELD"""
"""fix self interaction of water to value that reproduce the compressibility of pure water, u0 = 18.69kT, B = u0/(4 pi aev**2)**(3/2)"""
SysLoadFF = True
ForceFieldFile = 'ff_init.dat'

#fix Gaussian params of these pairs
fixPairs = ['HOH','HOH']
print('fixing Gaussian parameters of pairs {}'.format(fixPairs))

#Excluded volume size for each atom type: a_ev = 1/(number density of this CG atom type)
#aevs_self = {'Na+': 0.75/0.31, 'Cl-': 0.75/0.31, 'HOH': 0.75/0.31, 'A': 1.09/0.31,'A-': 1.09/0.31, 'B': 1.09/0.31, 'B+': 1.09/0.31, 'AE':1.09/0.31, 'AE-': 1.09/0.31, 'BE': 1.09/0.31, 'BE+':1.09/0.31}
aevs_self = {'Na+': 1., 'Cl-': 1., 'HOH': 1., 'Head': 0.5/0.31, 'Tail': 0.5/0.31}
for key in amino_acids:
    aevs_self.update({key: 0.5/0.31})
aCoul_self = aevs_self.copy()

#BondParams: (atom1,atom2):[Dist0,FConst,Label], FConts = kcal/mol/Angstrom**2
BondParams = {('All','All'):[0., 1.2, 'Bond']}
RLength_dict = {} #{'HOH': [r/lengthScale]}
#whether to fix a parameter [Dist0,FConst,Label]
IsFixedBond = {('All','All'):[True,False,True]}

#Pair interaction
#use spline or LJGauss for pair?
UseLJGauss = True
Cut = 10.0

#Gauss
#number of Gaussians for each pair type, if set {('All', 'All': n)}, all pairs have n Gaussians
NGaussDicts = {('All','All'): 1}

#Initial B
B0 = 0.5
u0_HOH_HOH = 15.1
BHOH_HOH = u0_HOH_HOH/(2 * np.pi * (aevs_self['HOH']**2+aevs_self['HOH']**2))**(3./2.) #B = u0/(2pi(ai^2+j^2))^3/2
LJGDist0 = 0.
LJGSigma = 1.
LJGEpsilon = 0.
#for now all pair types and all Gaussian interactions have the same fix
FixedB = False
FixedKappa = True

FixedDist0 = True
FixedSigma = True
FixedEpsilon = True

#Pair spline
PSplineParams = {}
PSplineNKnot = 10
NonbondEneSlopeInit = '1.kTperA'
FixedSpline = False

#Smeared Coul
SmearedCoulParams = {}
EwaldCoef = 2.4
#KMax = 10
SCoulShift = True
FixedCoef = False
FixedBornA = False
IsFixedCharge = True

#Ewald params
EwaldParams = {'ExcludeBondOrd': 0, 'Cut': Cut, 'Shift': True, 'Label': 'EW', 'Coef': EwaldCoef, 'EwaldNAtom': None, 'FixedCoef': FixedCoef} # 'KMax': KMax}

#External Potential, ExtPot can be list of many external potentials
UConst = 0.0 
ExtPot = {"UConst": UConst, "NPeriods": 1, "PlaneAxis": 2, "PlaneLoc": 0., 'AtomTypes':'A',"Label":"UExtSin"}

# Default Simulation Package Settings
import multiprocessing as mp
ncpu = mp.cpu_count()

sim.srel.base.n_process = min(6,ncpu)
sim.export.lammps.NeighOne = 8000
sim.export.lammps.UseTable2 = True
sim.export.lammps.InnerCutoff = 1.e-6
sim.export.lammps.NPairPotentialBins = 1000
sim.export.lammps.LammpsExec = 'lmp_omp'
sim.export.lammps.UseLangevin = True
sim.export.lammps.OMP_NumThread = 8 
sim.export.lammps.TableInterpolationStyle = 'spline' # More robust than spline for highly CG-ed systems
sim.srel.optimizetrajlammps.LammpsDelTempFiles = False
sim.srel.optimizetrajlammps.UseLangevin = True
sim.export.omm.platformName = 'OpenCL' # or 'OpenCL' or 'GPU' or 'CUDA'
sim.export.omm.device = -1 #-1 is default, let openmm choose its own platform.
sim.export.omm.NPairPotentialKnots = 500 #number of points used to spline-interpolate the potential
sim.export.omm.InnerCutoff = 0.001 #0.001 is default. Note that a small value is not necessary, like in the lammps export, because the omm export interpolates to zero
sim.srel.optimizetrajomm.OpenMMStepsMin = 0 #number of steps to minimize structure, 0 is default
sim.srel.optimizetrajomm.OpenMMDelTempFiles = False #False is Default
sim.export.omm.UseTabulated = True

"""---End of inputs---"""
        
"""Map AA trajectories to CG beads and get all unique CG atom types"""
#map AA trajectory if no CG trajs were provided
BoxLs = []
if len(CGtrajs) == 0:
    UniqueCGatomTypes = []
    for i, AAtraj in enumerate(AAtrajs):
        CGatomTypes, AAatomId, CGtraj, BoxL = mappoly.mapTraj(AAtraj,AAtops[i],nameMap, lengthScale, stride = strides[i])
        CGtrajs.append(CGtraj)
        BoxLs.append(BoxL)
        print "BoxL {}".format(BoxL)
        UniqueCGatomTypes.extend(CGatomTypes)
    UniqueCGatomTypes = np.unique(np.array(UniqueCGatomTypes)).tolist()
else:
    for i, CGtraj in enumerate(CGtrajs):
        CGtraj = pickleTraj(CGtraj)
        CGtrajs[i] = CGtraj
        BoxLs.append(CGtraj.FrameData['BoxL'])    
if UseLJGauss_CustomBead:
    _ = [x for y in ff_bead_map.values() for x in y]
    UniqueCGatomTypes = np.unique(_)
    print('UniqueCGatomTypes ',UniqueCGatomTypes)
"""Calculate forcefield parameters. This will likely to change"""
#get mixed term of aev and calculate kappa parameter for LJGauss, k = 1/(2 (ai^2+aj^2))
#Born radii = a_Coul * sqrt(pi)
aevs = {}
ks = {}
LJGaussParams = {} 
IsFixedLJGauss = {}
if UseLJGauss:
    #LJGauss param: (atom1,atom2):[B, kappa, Dist0, Cut, Sigma, Epsilon, Label ]
    if UseLJGauss_MolType:   #make pairs based on MolType, all beads on same MolType are treated the same
        for mol1 in MolNamesList[0]:
            for mol2 in MolNamesList[0]:
                atom1 = UniqueMolTypes[mol1][0] # pick first atom in this molecule to use for smearing length calculation 
                atom2 = UniqueMolTypes[mol2][0]
                a1 = aevs_self[atom1]
                a2 = aevs_self[atom2]
                a12 =  np.sqrt((a1**2 + a2**2)/2)
                kappa = 1/(2*(a1**2 + a2**2))
                aevs.update({(atom1,atom2): a12})
                ks.update({(atom1,atom2): kappa})
                pair = ('MOL_{}'.format(mol1), 'MOL_{}'.format(mol2))
                if not pair in  LJGaussParams.keys() and not (pair[1],pair[0]) in LJGaussParams.keys():
                    if pair == ('MOL_HOH','MOL_HOH'):
                        print('fix params of pair {}'.format(pair))
                        LJGaussParams.update({pair: [BHOH_HOH, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(mol1,mol2)]})
                        IsFixedLJGauss.update({pair: [True, True, FixedDist0, FixedSigma, FixedEpsilon, True]})
                    elif (mol1,mol2) in fixPairs or (mol2,mol1) in fixPairs:
                        print('fix params of pair {}'.format(pair))
                        LJGaussParams.update({pair: [B0, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(mol1,mol2)]})
                        IsFixedLJGauss.update({pair: [True, True, True, True, True, True]})
                    else:
                        LJGaussParams.update({pair: [B0, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(mol1,mol2)]})
                        IsFixedLJGauss.update({pair: [FixedB, FixedKappa, FixedDist0, FixedSigma, FixedEpsilon, True]})
    elif UseLJGauss_CustomBead:
        LJGauss_groups = [x for x in ff_bead_map.keys()]
        for i in range(len(LJGauss_groups)):
            for j in range(len(LJGauss_groups)):
                group1 = LJGauss_groups[i]
                group2 = LJGauss_groups[j]
                atoms1 = tuple(ff_bead_map[group1])
                atoms2 = tuple(ff_bead_map[group2])
                # take average radius if atomsi include many bead types
                a1 = np.mean([aevs_self[atom1] for atom1 in atoms1])
                a2 = np.mean([aevs_self[atom2] for atom2 in atoms2]) 
                a12 =  np.sqrt((a1**2 + a2**2)/2)
                kappa = 1/(2*(a1**2 + a2**2))
                aevs.update({(group1,group2): a12})
                ks.update({(group1,group2): kappa})
                #add param if this pair has not been added to param dictionary
                pot_name = 'LJGauss{}_{}'.format(group1,group2)
                if (atoms1,atoms2) in LJGaussParams.keys() or (atoms2,atoms1) in LJGaussParams.keys(): continue
                if (group1,group2) == ('HOH','HOH'):
                    print('fix params of pair {}'.format((group1,group2)))
                    LJGaussParams.update({(atoms1,atoms2): [BHOH_HOH, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, pot_name]})
                    IsFixedLJGauss.update({(atoms1,atoms2): [True, True, FixedDist0, FixedSigma, FixedEpsilon, True]})
                elif (group1,group2) in fixPairs or (group2,group1) in fixPairs:
                    print('fix params of pair {}'.format((group1,group2)))
                    LJGaussParams.update({(atoms1,atoms2): [B0, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, pot_name]})
                    IsFixedLJGauss.update({(atoms1,atoms2): [True, True, True, True, True, True]})

                else:
                    LJGaussParams.update({(atoms1,atoms2): [B0, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, pot_name]})
                    IsFixedLJGauss.update({(atoms1,atoms2): [FixedB, FixedKappa, FixedDist0, FixedSigma, FixedEpsilon, True]})
        print('IsFixedLJGauss ',IsFixedLJGauss)
    else: # make pairs based on AtomType
        for i in range(len(UniqueCGatomTypes)):
            for j in range(len(UniqueCGatomTypes)):
                atom1 = UniqueCGatomTypes[i]
                atom2 = UniqueCGatomTypes[j]
                a1 = aevs_self[atom1]
                a2 = aevs_self[atom2]
                a12 =  np.sqrt((a1**2 + a2**2)/2)
                kappa = 1/(2*(a1**2 + a2**2))
                aevs.update({(atom1,atom2): a12})
                ks.update({(atom1,atom2): kappa})
                #add param if this pair has not been added to param dictionary
                if not (atom1,atom2) in LJGaussParams.keys() and not (atom2,atom1) in LJGaussParams.keys():
                    if (atom1,atom2) == ('HOH','HOH'):
                        print('fix params of pair {}'.format((atom1,atom2)))
                        LJGaussParams.update({(atom1,atom2): [BHOH_HOH, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(atom1,atom2)]})
                        IsFixedLJGauss.update({(atom1,atom2): [True, True, FixedDist0, FixedSigma, FixedEpsilon, True]})
                    elif (atom1,atom2) in fixPairs or (atom2,atom1) in fixPairs:
                        print('fix params of pair {}'.format((atom1,atom2)))
                        LJGaussParams.update({(atom1,atom2): [B0, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(atom1,atom2)]})
                        IsFixedLJGauss.update({(atom1,atom2): [True, True, True, True, True, True]}) 
                    else:
                        LJGaussParams.update({(atom1,atom2): [B0, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(atom1,atom2)]})
                        IsFixedLJGauss.update({(atom1,atom2): [FixedB, FixedKappa, FixedDist0, FixedSigma, FixedEpsilon, True]})
    f = open('aevs.txt','w')
    f.write('{}'.format(aevs))
    f.close()
else:
    #PSpline param : (atom1, atom2): [NKnot, Cut, NonbondEneSlopeInit,  Label]
    for i in range(len(UniqueCGatomTypes)):
        for j in range(len(UniqueCGatomTypes)):
            atom1 = UniqueCGatomTypes[i]
            atom2 = UniqueCGatomTypes[j]
            if not (atom1,atom2) in PSplineParams.keys() and not (atom2,atom1) in PSplineParams.keys():
                PSplineParams.update({(atom1,atom2): [PSplineNKnot, Cut, NonbondEneSlopeInit, 'PSpline{}_{}'.format(atom1,atom2), FixedSpline]})
#SmearedCoulParams: (atom1,atom2): [BornA, Cut, Shift, FixedCoef,FixedBornA, Label, Coef]
BornAs = {}
for i in range(len(UniqueCGatomTypes)):
    for j in range(len(UniqueCGatomTypes)):
        atom1 = UniqueCGatomTypes[i]
        atom2 = UniqueCGatomTypes[j]
        a1 = aCoul_self[atom1]
        a2 = aCoul_self[atom2]
        a12 =  np.sqrt((a1**2 + a2**2)/2)
        BornA = a12*np.sqrt(np.pi)
        BornAs.update({(atom1,atom2): BornA})
        if all([charges[atom1], charges[atom2]]) != 0.:
            if not (atom1,atom2) in SmearedCoulParams.keys() and not (atom2,atom1) in SmearedCoulParams.keys():
                SmearedCoulParams.update({(atom1,atom2):[BornA, Cut, SCoulShift, FixedCoef, FixedBornA, 'SmearCoul{}_{}'.format(atom1,atom2), EwaldCoef]})    

"""Create World to share among Systems"""
World = system.CreateWorld(UniqueCGatomTypes, UniqueMolTypes, charges, elements, Units=Units, RLength_dict=RLength_dict)
if SrelOnNeutralSys:
    WorldNeutral = system.CreateWorld(UniqueCGatomTypes, UniqueMolTypes, chargesNeutral, elements, Units=Units, RLength_dict=RLength_dict)
"""Create systems and optimizers"""
Systems = []
Maps = []
Opts = []
NAtoms = []
NMols = []
if len(PresSet) > 0:
        UseNPT = True
        IntParams.update({'ensemble':'NPT'})
else:
        UseNPT = False
        IntParams.update({'ensemble':'NVT'})
if SysLoadFF:
    ForceFieldFile = ForceFieldFile
else:
    ForceFieldFile = None
for i, CGtraj in enumerate(CGtrajs):
    NMolsDict = NMolsDicts[i]
    MolNames = MolNamesList[i]
    Temp = TempSet[i]
    SysName = Name+str(i)
    BoxL = BoxLs[i] 
    MolIdRgs = MolIdRgList[i]
    RgTars = RgTarList[i]
    if UseNPT:
        Pres = PresSet[i]
    else:
        Pres = None
    if len(StepScales) > 0:
        StepScale = StepScales[i]
    else:
        StepScale = None
    #create system and add forcefield, then create optimizer for each system
    if SrelOnNeutralSys:
        print('Optimizing params in neutral system but run MD on system with full electrostatics')
        ElecSys,_ = system.CreateSystem(SysName+'Elec', World, BoxL, MolNames, NMolsDict, IsFixedCharge, Temp, Pres, IntParams,ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot,
                              RandomMul=RandomMul)
        Sys,measureRgs  = system.CreateSystem(SysName+'Neutral', WorldNeutral, BoxL, MolNames, NMolsDict, IsFixedCharge, Temp, Pres, IntParams,ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot,
                              RandomMul=RandomMul)
    else:
        Sys,measureRgs = system.CreateSystem(SysName, World, BoxL, MolNames, NMolsDict, IsFixedCharge, Temp, Pres, IntParams,ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot,
                              RandomMul=RandomMul)
        ElecSys = None

    ParamString = Sys.ForceField.ParamString()
    paramfile = open(SysName+'_ff.dat','w')
    paramfile.write(ParamString)
    paramfile.close()

    if ElecSys:
        ParamString = ElecSys.ForceField.ParamString()
        paramfile = open(SysName+'_ele_ff.dat','w')
        paramfile.write(ParamString)
        paramfile.close()

    Opt = optimizer.CreateOptimizer(Sys, CGtraj, UseLammps, UseOMM, UseSim, StepsEquil, StepsProd, StepsStride, StepScale, UseWPenalty,ElecSys = ElecSys, RgConstrain = RgConstrain,RgTars=RgTars,measureRgs=measureRgs, ConstrainMode=ConstrainMode)

    Opts.append(Opt)
    Systems.append(Sys)
    NAtoms.append(Sys.NAtom)
    NMols.append(Sys.NMol)

"""Run Optimzation"""
# Just always using the OptimizeMultiTrajClass
Weights = [1.]*len(Opts)

if WeightSysByNMol:
    Weights = []
    for NMol in NMols:
        Weights.append(np.max(NMols)/float(NMol))
    print ('Weights for Expanded Ensemble are:')
    print (Weights)
    
elif WeightSysByNAtom:
    Weights = []
    for NAtom in NAtoms:
        Weights.append(np.max(NAtoms)/float(NAtom))
    print ('Weights for Expanded Ensemble are:')
    print (Weights)
optimizer.RunOpt(Opts, Weights, Name, UseWPenalty, MaxIter, SteepestIter, RgConstrain, StageCoefs=StageCoefs)

