#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:11:03 2019

@author: my
"""
import sim, pickleTraj
import mappoly, system, optimizer
import numpy as np
"""assume all atom indices in trajectory follow the same order as definition of MolTypes, NMons, and NMols
To do:
   allow mapping other than 1:1 for polymer
   read bond info with mdtraj to calculate number of monomers (NMons) ber bead and NMols
   check read BoxL from AAtraj, unit consistency
   make sure all systems in the EE have same forcefield
   how to run neutral system along with charged system? lammps doesn't allow coulomb for neutral beads
   lammps: fix overlay pair style coeff, use table2 when appropriate"""

"""---Inputs---"""  
kB = 1.380658e-23
TempRef = 298.15
kT = kB*TempRef
kTkJmol = kT/1000.0*6.022e23
kTkcalmol = kTkJmol/4.184

"""TOPOLOGY"""
AAtrajs = ['trajectory_xp0.1_N12_f0_V157_LJPME_298K_NVT_Uext0.dcd']
AAtops = ['AA12_f0_opc_gaff2_w0.13.parm7']
stride = 10
#map AA residue to CG bead name
nameMap = {'Na+':'Na+', 'Cl-':'Cl-', 'HOH': 'HOH', 'WAT': 'HOH',
               'ATP':'A', 'AHP':'A', 'AP': 'A', 'ATD': 'A-', 'AHD': 'A-', 'AD': 'A-',
               'NTP':'B+', 'NHP':'B+', 'NP': 'B+', 'NTD': 'B', 'NHD': 'B', 'ND': 'B'}
CGtrajs = []
#CGtrajs = ['trajectory_xp01_N12_f0_V157_LJPME_298K_NVT_Uext0_mapped.lammpstrj.gz']
#provide UniqueCGatomTypes if CGtrajs is not an empty list
UniqueCGatomTypes = ['A','HOH']

#name of molecules in systems
#must in the right sequence as molecules in the trajectory
MolNamesList = [['PAA','HOH']]
# nSys x molecule types   
MolTypesDicts = [{'PAA':['A']*12,'Na+':['Na+'],'Cl-':['Cl-'],'HOH':['HOH']}]
# number of molecules for each molecule type, nSys x molecule types
NMolsDicts = [{'PAA':15,'Na+':0,'Cl-':0,'HOH':4728}]
charges = {'Na+': 1., 'Cl-': -1., 'HOH': 0., 'A': 0,'A-': -1., 'B': 0., 'B-': 1.}

Name = 'PAA'

"""INTEGRATION PARAMS"""
dt = 0.005 #ps
TempSet = [298.15] #K
PresSet = [] #enter values to enable NPT

PresSet = np.array(PresSet)
PresSet *= 101325. #joules/m**3
PresSet *= (10.**-10.)**3. * 0.000239006 *6.022e23 #kcal/mol/A**3

IntParams = {'TimeStep': dt, 'LangevinGamma': 1/(100*dt)}

"""SREL OPT"""

UseLammps = True
UseOMM = False
UseSim = False
ScaleRuns = True
StepScales = [] #set to empty if don't want to scale steps
StepsEquil = 10000
StepsProd = 200000
StepsStride = 100
WeightSysByNMol = False
WeightSysByNAtom = True

UseWPenalty = False
StageCoefs = []

MaxIter=None
SteepestIter=0

"""FORCEFIELD"""
"""fix self interaction of water to value that reproduce the compressibility of pure water, u0 = 18.69kT, B = u0/(4 pi aev**2)**(3/2)"""
SysLoadFF = False
ForceFieldFile = 'ff.dat'
#Excluded volume size for each atom type: a_ev = 1/(number density of this CG atom type)
aevs_self = {'Na+': 1.8, 'Cl-': 1.8, 'HOH': 3.1, 'A': 4.5,'A-': 4.5, 'B': 4.5, 'B-': 4.5} #Angstrom
aCoul_self = aevs_self.copy()
#BondParams: (atom1,atom2):[Dist0,FConst,Label], FConts = kcal/mol/Angstrom**2
BondParams = {('A','A'):[4., 100*kTkcalmol, 'BondA_A']}
#{('A','A-'):[4., 50*kTkcalmol, 'BondA_A-'], ('A','A'):[4., 50*kTkcalmol, 'BondA_A'],
#               ('B','B+'):[4., 50*kTkcalmol, 'BondB_B+'], ('B','B'):[4., 50*kTkcalmol, 'BondB_B']}
#whether to fix a parameter
IsFixedBond = {('A','A-'):[False,False,True], ('A','A'):[False,False,True], ('A-','A-'):[False,False,True],
               ('B','B+'):[False,False,True], ('B','B'):[False,False,True], ('B+','B+'):[False,False,True]}
#set cut to be 5 * the largest aev
Cut = 5 * np.max(aevs_self.values())
#Initial B
B0 = 20. * kTkcalmol
#BHOH_HOH = 18.69 * kTkcalmol
LJGDist0 = 0.
LJGSigma = 1.
LJGEpsilon = 0.
LJGaussParams = {} 
IsFixedLJGauss = {}
#for now all pairs have same fixed paramters
FixedB = False
FixedKappa = True
FixedDist0 = True
FixedSigma = True
FixedEpsilon = True
#Smeared Coul
SmearedCoulParams = {}
SCoulShift = True
FixedCoef = False
FixedBornA = True
IsFixedCharge = True
#Ewald params
EwaldParams = {'ExcludeBondOrd': 0, 'Cut': Cut, 'Shift': True, 'Label': 'EW'}

#External Potential
UConst = 0.*kTkcalmol
ExtPot = {"UConst": UConst, "NPeriods": 1, "PlaneAxis": 2, "PlaneLoc": 0., 'AtomTypes':['HOH']}

# Default Simulation Package Settings
sim.export.lammps.UseTable2 = True
sim.export.lammps.InnerCutoff = 1.e-6
sim.export.lammps.NPairPotentialBins = 1000
sim.export.lammps.LammpsExec = 'lmp_omp'
sim.export.lammps.UseLangevin = True
sim.export.lammps.OMP_NumThread = 5
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
        CGatomTypes, AAatomId, CGtraj, BoxL = mappoly.mapTraj(AAtraj,AAtops[i],nameMap, stride = stride)
        CGtrajs.append(CGtraj)
        BoxLs.append(BoxL)
        print "BoxL {}".format(BoxL)
        UniqueCGatomTypes.append(CGatomTypes)
    UniqueCGatomTypes = np.unique(np.array(UniqueCGatomTypes))
else:
    for i, CGtraj in enumerate(CGtrajs):
        CGtraj = pickleTraj(CGtraj)
        CGtrajs[i] = CGtraj
        BoxLs.append(CGtraj.FrameData['BoxL'])    
"""Calculate forcefield parameters. This will likely to change"""
#get mixed term of aev and calculate kappa parameter for LJGauss, k = 1/(4aev^2)
#Born radii = a_Coul * sqrt(pi)
aevs = {}
ks = {}
#LJGauss param: (atom1,atom2):[B, kappa, Dist0, Cut, Sigma, Epsilon, Label ]
for i in range(len(UniqueCGatomTypes)):
    for j in range(len(UniqueCGatomTypes)):
        atom1 = UniqueCGatomTypes[i]
        atom2 = UniqueCGatomTypes[j]
        a1 = aevs_self[atom1]
        a2 = aevs_self[atom2]
        amean = np.mean([a1,a2])
        kappa = 1/(4*amean**2)
        aevs.update({(atom1,atom2): amean})
        ks.update({(atom1,atom2): 1/(4*amean**2)})
        #add param if this pair has not been added to param dictionary
        if not (atom1,atom2) in LJGaussParams.keys() and not (atom2,atom1) in LJGaussParams.keys():
#            if (atom1,atom2) == ('HOH','HOH'):
#                LJGaussParams.update({(atom1,atom2): [BHOH_HOH, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(atom1,atom2)]})
#                IsFixedLJGauss.update({(atom1,atom2): [True, True, FixedDist0, FixedSigma, FixedEpsilon, True]})
#            else:
             LJGaussParams.update({(atom1,atom2): [B0, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(atom1,atom2)]})
             IsFixedLJGauss.update({(atom1,atom2): [FixedB, FixedKappa, FixedDist0, FixedSigma, FixedEpsilon, True]})
f = open('aevs.txt','w')
f.write('{}'.format(aevs))
f.close()

#SmearedCoulParams: (atom1,atom2): [BornA, Cut, Shift, FixedCoef,FixedBornA, Label]
BornAs = {}
for i in range(len(UniqueCGatomTypes)):
    for j in range(len(UniqueCGatomTypes)):
        atom1 = UniqueCGatomTypes[i]
        atom2 = UniqueCGatomTypes[j]
        a1 = aCoul_self[atom1]
        a2 = aCoul_self[atom2]
        amean = np.mean([a1,a2])
        BornA = amean*np.sqrt(np.pi)
        BornAs.update({(atom1,atom2): BornA})
        if not (atom1,atom2) in SmearedCoulParams.keys() and not (atom2,atom1) in SmearedCoulParams.keys():
            SmearedCoulParams.update({(atom1,atom2):[BornA, Cut, SCoulShift, FixedCoef, FixedBornA, 'SmearCoul{}_{}'.format(atom1,atom2)]})    

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
for i, MolTypesDict in enumerate(MolTypesDicts):
    NMolsDict = NMolsDicts[i]
    MolNames = MolNamesList[i]
    Temp = TempSet[i]
    SysName = Name+str(i)
    CGtraj = CGtrajs[i]
    BoxL = BoxLs[i] 
    if UseNPT:
        Pres = PresSet[i]
    else:
        Pres = None
    if len(StepScales) > 0:
        StepScale = StepScales[i]
    else:
        StepScale = None
    #create system and add forcefield, then create optimizer for each system
    Sys = system.CreateSystem(SysName, BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams,ForceFieldFile,
                              LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, ExtPot)
    Opt = optimizer.CreateOptimizer(Sys, CGtraj, UseLammps, UseOMM, UseSim, StepsEquil, StepsProd, StepsStride, StepScale, UseWPenalty)

    Opts.append(Opt)
    Systems.append(Sys)
    NAtoms.append(Sys.NAtom)
    NMols.append(Sys.NMol)

"""Run Optimization"""
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

Optimizer = sim.srel.OptimizeMultiTrajClass(Opts, Weights=Weights)
Optimizer.FilePrefix = ("{}".format(Name))

if not UseWPenalty:
    Optimizer.RunConjugateGradient(MaxIter=MaxIter, SteepestIter=SteepestIter)
else:
    Optimizer.RunStages(StageCoefs = StageCoefs)
