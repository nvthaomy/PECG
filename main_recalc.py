#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:11:03 2019

@author: my
"""
import sys
sys.path.append('/home/mnguyen/bin/PECG/')

import sim, pickleTraj
#import mappoly_implicit as mappoly
import mappoly
import system, optimizer
import numpy as np
import argparse as ap

parser = ap.ArgumentParser(description='recalculate hessian')
parser.add_argument('ff_file',type=str, default=None, help='force field parameter filename')
#parser.add_argument('label',type=str, help='label to distinguish output from this calculation')
#parser.add_argument('-test', action='store_true', help='whether to use shorter test trajectory')
#parser.add_argument('-recalc',type=bool, default=True, help='whether or not we are just recalculating')
args = parser.parse_args()


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
lengthScale  = 3.12 # = CG length scale/ AA length scale = 3.1 A/ 1 A
if Units == sim.units.DimensionlessUnits:
    lengthScale = lengthScale
else:
    lengthScale = 1.
#dimensionless scales:
# length: a_water = 3.1 A
# energy: kT
# pressure: kT/a_water**3

"""TOPOLOGY"""
#map AA traj to CG traj
AAtrajs = ['/home/mnguyen/PE/nacl/AA/2Mnacl_opc/NVT/uext2_nacl/trajectory_2Mnacl_opc_Uext2_nacl_298K.dcd']
AAtops = ['/home/mnguyen/PE/nacl/AA/2Mnacl_opc/NVT/uext2_nacl/6432opc_245nacl.parm7']
stride = 10

#map AA residue to CG bead name
nameMap = {'Na+':'Na+', 'Cl-':'Cl-', 'HOH': 'HOH', 'WAT': 'HOH',
               'ATP':'A', 'AHP':'A', 'AP': 'A', 'ATD': 'A-', 'AHD': 'A-', 'AD': 'A-',
               'NTP':'B+', 'NHP':'B+', 'NP': 'B+', 'NTD': 'B', 'NHD': 'B', 'ND': 'B'}

#if use different bead type for end monomers:
#nameMap = {'Na+':'Na+', 'Cl-':'Cl-', 'HOH': 'HOH', 'WAT': 'HOH',
#               'ATP':'AE', 'AHP':'AE', 'AP': 'A', 'ATD': 'AE-', 'AHD': 'AE-', 'AD': 'A-',
#               'NTP':'BE+', 'NHP':'BE+', 'NP': 'B+', 'NTD': 'BE', 'NHD': 'BE', 'ND': 'B'}

#CG trajectories, if CGtrajs = [] will map AA traj, else will use provided CGtrajs
CGtrajs = []
#provide UniqueCGatomTypes if CGtrajs is not empty
UniqueCGatomTypes = ['Na+','Cl-','HOH']

#name of molecules in systems
#must in the right sequence as molecules in the trajectory
MolNamesList = [['Na+','Cl-','HOH']]

# Topology of molecules in all systems: nSys x molecule types   
MolTypesDicts = [{'PAA':['A','A']*6,'Na+':['Na+'],'Cl-':['Cl-'],'HOH':['HOH']}]

# number of molecules for each molecule type, nSys x molecule types
NMolsDicts = [{'PAA':0,'Na+': 245,'Cl-': 245,'HOH': 6432}]

charges         = {'Na+': 1., 'Cl-': -1., 'HOH': 0., 'A': 0, 'A-': -1., 'B': 0., 'B+': 1., 'AE':0., 'BE':0., 'AE-': -1., 'BE+':1.}
chargesNeutral  = {'Na+': 0., 'Cl-': 0.,  'HOH': 0., 'A': 0, 'A-': 0.,  'B': 0., 'B+': 0., 'AE':0., 'BE':0., 'AE-': 0. , 'BE+':0.}

Name = 'nacl'

"""INTEGRATION PARAMS"""
#real units: dt (ps), temp (K), pressure (atm)
dt = 0.1
#TempSet and PresSet are lists 
TempSet = [1.]
PresSet = [] #enter values to enable NPT

PresSet = np.array(PresSet)
#PresSet *= 101325. #joules/m**3
#PresSet *= (10.**-10.)**3. * 0.000239006 *6.022e23 #kcal/mol/A**3

IntParams = {'TimeStep': dt, 'LangevinGamma': 1/(100*dt)}

"""SREL OPT"""
SrelOnNeutralSys = True
UseLammps = False
UseOMM = True
UseSim = False
ScaleRuns = True
StepScales = [] #set to empty if don't want to scale steps
StepsEquil = 200./dt
StepsProd = 10000./dt
StepsStride = 100
WeightSysByNMol = False
WeightSysByNAtom = True

UseWPenalty = False
StageCoefs = []

MaxIter=None
SteepestIter=0

"""FORCEFIELD"""
"""fix self interaction of water to value that reproduce the compressibility of pure water, u0 = 18.69kT, B = u0/(4 pi aev**2)**(3/2)"""
SysLoadFF = True
ForceFieldFile = args.ff_file

#fix Gaussian params of these pairs
fixPairs = [('HOH','HOH')]
print('fixing Gaussian parameters of pairs {}'.format(fixPairs))

#Excluded volume size for each atom type: a_ev = 1/(number density of this CG atom type)
aevs_self = {'Na+': 1., 'Cl-': 1., 'HOH': 1., 'A': 4.5/3.12,'A-': 4.5/3.12, 'B': 4.5/3.12, 'B+': 4.5/3.12, 'AE':4.5/3.12, 'AE-': 4.5/3.12, 'BE': 4.5/3.12, 'BE+':4.5/3.12}
aCoul_self = aevs_self.copy()

#BondParams: (atom1,atom2):[Dist0,FConst,Label], FConts = kcal/mol/Angstrom**2
BondParams = {} #{('A-','A-'):[1., 2000, 'BondA-_A-']}
#{('A','A-'):[4., 50*kTkcalmol, 'BondA_A-'], ('A','A'):[4., 50*kTkcalmol, 'BondA_A'],
#               ('B','B+'):[4., 50*kTkcalmol, 'BondB_B+'], ('B','B'):[4., 50*kTkcalmol, 'BondB_B']}

#whether to fix a parameter [Dist0,FConst,Label]
IsFixedBond = {('A','A-'):[False,False,True], ('A','A'):[False,False,True], ('A-','A-'):[False,False,True],
               ('B','B+'):[False,False,True], ('B','B'):[False,False,True], ('B+','B+'):[False,False,True],
               ('A','AE'):[False,False,True], ('A','AE-'):[False,False,True], ('A-','AE'):[False,False,True], ('A-','AE-'):[False,False,True],
               ('B','BE'):[False,False,True], ('B','BE+'):[False,False,True], ('B+','BE'):[False,False,True], ('B+','BE+'):[False,False,True]}
#Pair interaction
#use spline or LJGauss for pair?
UseLJGauss = True

Cut = 6.

#Gauss
#number of Gaussians for each pair type, if set {('All', 'All': n)}, all pairs have n Gaussians
NGaussDicts = {('All','All'): 1}

#Initial B
B0 = 0.5
u0_HOH_HOH = 15.1
BHOH_HOH = 0.1 #u0_HOH_HOH/(2 * np.pi * (aevs_self['HOH']**2+aevs_self['HOH']**2))**(3./2.) #B = u0/(2pi(ai^2+j^2))^3/2
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
SCoulShift = True
FixedCoef = False
FixedBornA = False
IsFixedCharge = True

#Ewald params
EwaldParams = {'ExcludeBondOrd': 0, 'Cut': Cut, 'Shift': True, 'Label': 'EW', 'Coef': EwaldCoef, 'EwaldNAtom': None, 'FixedCoef': FixedCoef}

#External Potential, ExtPot can be list of many external potentials
UConst = 2.0 
ExtPot0 = {"UConst": UConst, "NPeriods": 1, "PlaneAxis": 2, "PlaneLoc": 0., 'AtomTypes':'Na+',"Label":"UExtSin_Na"} 
ExtPot1 = {"UConst": UConst, "NPeriods": 1, "PlaneAxis": 2, "PlaneLoc": 0., 'AtomTypes':'Cl-',"Label":"UExtSin_Cl"} 
ExtPot = [ExtPot0,ExtPot1]
IsFixedExtPot =  {"UConst": False, "NPeriods":True}
 
# Default Simulation Package Settings
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
        CGatomTypes, AAatomId, CGtraj, BoxL = mappoly.mapTraj(AAtraj,AAtops[i],nameMap, lengthScale, stride = stride)
        CGtrajs.append(CGtraj)
        BoxLs.append(BoxL)
        print "BoxL {}".format(BoxL)
        UniqueCGatomTypes.extend(CGatomTypes)
    UniqueCGatomTypes = np.unique(np.array(UniqueCGatomTypes))
else:
    for i, CGtraj in enumerate(CGtrajs):
        CGtraj = pickleTraj(CGtraj)
        CGtrajs[i] = CGtraj
        BoxLs.append(CGtraj.FrameData['BoxL'])    
"""Calculate forcefield parameters. This will likely to change"""
#get mixed term of aev and calculate kappa parameter for LJGauss, k = 1/(2 (ai^2+aj^2))
#Born radii = a_Coul * sqrt(pi)
aevs = {}
ks = {}
LJGaussParams = {} 
IsFixedLJGauss = {}
if UseLJGauss:
    #LJGauss param: (atom1,atom2):[B, kappa, Dist0, Cut, Sigma, Epsilon, Label ]
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
    if SrelOnNeutralSys:
        print('Optimizing params in neutral system but run MD on system with full electrostatics')
        ElecSys, _ = system.CreateSystem(SysName+'Elec', BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams,ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, IsFixedExtPot=IsFixedExtPot,Units = Units)
        Sys, measureRgs = system.CreateSystem(SysName+'Neutral', BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, chargesNeutral, IsFixedCharge, Temp, Pres, IntParams,ForceFieldFile,          
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, IsFixedExtPot=IsFixedExtPot, Units = Units) 
    else:
        Sys, measureRgs = system.CreateSystem(SysName, BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams,ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, IsFixedExtPot=IsFixedExtPot, Units = Units)
        ElecSys = None

    Opt = optimizer.CreateOptimizer(Sys, CGtraj, UseLammps, UseOMM, UseSim, StepsEquil, StepsProd, StepsStride, StepScale, UseWPenalty,ElecSys = ElecSys, recalc=True)

    Opts.append(Opt)
    Systems.append(Sys)
    NAtoms.append(Sys.NAtom)
    NMols.append(Sys.NMol)

"""Run Optimzation"""
optimizer.recalc(Opts,Prefix=Name+'_recalc')

