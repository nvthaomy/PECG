#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:11:03 2019

@author: my
"""
import sim, pickleTraj
import sys
sys.path.append('/home/mnguyen/bin/PECG/')
import system, optimizer, fepcg
import numpy as np
import MDAnalysis as mda
import os,sys
cwd = os.getcwd()

"""assume all atom indices in trajectory follow the same order as definition of MolTypes, NMons, and NMols
To do:
   allow mapping other than 1:1 for polymer
   read bond info with mdtraj to calculate number of monomers (NMons) ber bead and NMols
   check read BoxL from AAtraj, unit consistency
   make sure all systems in the EE have same forcefield
   how to run neutral system along with charged system? lammps doesn't allow coulomb for neutral beads
   lammps: fix overlay pair style coeff, use table2 when appropriate
   
   generalize to calculate Rg if there are two different polymers, currently assume calculate Rg on PAA"""

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
nameMap = {'Na+':'Na+', 'Cl-':'Cl-', 'HOH': 'HOH', 'WAT': 'HOH',
               'ATP':'A', 'AHP':'A', 'AP': 'A', 'ATD': 'A-', 'AHD': 'A-', 'AD': 'A-',
               'NTP':'B+', 'NHP':'B+', 'NP': 'B+', 'NTD': 'B', 'NHD': 'B', 'ND': 'B'}
#dimensionless (average volume from NPT)
L = __vol__**(1./3.) 
BoxLs = [[L,L,L]]

#name of molecules in systems
#must in the right sequence as molecules in the trajectory
MolNamesList = [__MolNamesList__]
# nSys x molecule types   
MolTypesDicts = [{'PAA':__PAAstructure__,'Na+':['Na+'],'Cl-':['Cl-'],'HOH':['HOH']}]
# number of molecules for each molecule type, nSys x molecule types
NMolsDicts = [{'PAA':__nPAA__,'Na+':__nNa__,'Cl-':__nCl__,'HOH':__nHOH__}]
charges         = {'Na+': 1., 'Cl-': -1., 'HOH': 0., 'A': 0, 'A-': -1., 'B': 0., 'B+': 1., 'AE':0., 'BE':0., 'AE-': -1., 'BE+':1.}

Name = __Name__

#get UniqueCGatomTypes
UniqueCGatomTypes = []
for i, NMolsDict in enumerate(NMolsDicts):
    MolTypesDict = MolTypesDicts[i]
    for Mol,NMol in NMolsDict.items():
        UniqueCGatomTypes.extend(NMol*MolTypesDict[Mol]) 
UniqueCGatomTypes = np.unique(np.array(UniqueCGatomTypes))
print(UniqueCGatomTypes) 

"""INTEGRATION PARAMS"""
#real units: dt (ps), temp (K), pressure (atm)
dt = __dt__
TempSet = [1.]
PresSet = [__P__] #enter values to enable NPT

PresSet = np.array(PresSet)
#PresSet *= 101325. #joules/m**3
#PresSet *= (10.**-10.)**3. * 0.000239006 *6.022e23 #kcal/mol/A**3

IntParams = {'TimeStep': dt, 'LangevinGamma': 1/(100*dt)}

"""MD OPT"""
RunMD = True
UseLammps = __UseLammps__
UseOMM = __UseOMM__
UseSim = False
StepsMin = 1000
StepScales = [] #set to empty if don't want to scale steps
StepsEquil = __equilTau__/dt 
StepsProd = __tau__/dt
StepsStride = __Stride__
MDRestartFile = None #None: dont read restart file

"""FEP Params"""
CalChemPot = False
FEPDir = __FEPDir__
FEPMolNames = __FEPMolNames__ #['HOH']
ThermoSlice = __ThermoSlice__ 
TrajSlice = __TrajSlice__
Warmup = 50 
nDraw = 1
nInsert = __nInsert__
nDelete = __nDelete__
SetReRunRefState = False #always rerun traj if UseLammps to make sure potential energy is in the correct unit

"""FORCEFIELD"""
"""fix self interaction of water to value that reproduce the compressibility of pure water, u0 = 18.69kT, B = u0/(4 pi aev**2)**(3/2)"""
SysLoadFF = True #always read force field when run MD
ForceFieldFile = __ff__

#Excluded volume size for each atom type: a_ev = 1/(number density of this CG atom type)
aevs_self = {'Na+': 1., 'Cl-': 1., 'HOH': 1., 'A': 4.5/3.1,'A-': 4.5/3.1, 'B': 4.5/3.1, 'B+': 4.5/3.1, 'AE':4.5/3.1, 'AE-': 4.5/3.1, 'BE': 4.5/3.1, 'BE+':4.5/3.1}
aCoul_self = aevs_self.copy()

#BondParams: (atom1,atom2):[Dist0,FConst,Label], FConts = kcal/mol/Angstrom**2
BondParams = {}
#               ('B','B+'):[4., 50*kTkcalmol, 'BondB_B+'], ('B','B'):[4., 50*kTkcalmol, 'BondB_B']}
#whether to fix a parameter
IsFixedBond = {('A','A-'):[False,False,True], ('A','A'):[False,False,True], ('A-','A-'):[False,False,True],
               ('B','B+'):[False,False,True], ('B','B'):[False,False,True], ('B+','B+'):[False,False,True],
               ('A','AE'):[False,False,True], ('A','AE-'):[False,False,True], ('A-','AE'):[False,False,True], ('A-','AE-'):[False,False,True],
               ('B','BE'):[False,False,True], ('B','BE+'):[False,False,True], ('B+','BE'):[False,False,True], ('B+','BE+'):[False,False,True]}

#Pair interaction
#use spline or LJGauss for pair?
UseLJGauss = True

#set cut to be 5 * the largest aev
Cut = __cut__ #5 * np.max(aevs_self.values())

#Gauss
#number of Gaussians for each pair type
NGaussDicts = {('All','All'): 1}

#Initial B
B0 = 0.5
u0_HOH_HOH = 15.1
BHOH_HOH = u0_HOH_HOH/(2 * np.pi * (aevs_self['HOH']**2+aevs_self['HOH']**2))**(3./2.) #B = u0/(2pi(ai^2+j^2))^3/2

LJGDist0 = 0.
LJGSigma = 1.
LJGEpsilon = 0.
LJGaussParams = {} 
IsFixedLJGauss = {}
#for now all pair types and all Gaussian interactions have same fixed paramters
FixedB = False
FixedKappa = False

FixedDist0 = True
FixedSigma = True
FixedEpsilon = True

#Pair spline
PSplineParams = {}
PSplineNKnot = 20
NonbondEneSlopeInit = '1.kTperA'

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
UConst = 0.
ExtPot = {"UConst": UConst, "NPeriods": 1, "PlaneAxis": 2, "PlaneLoc": 0., 'AtomTypes':'HOH',"Label":"UExtSin"}

# Default Simulation Package Settings
sim.export.lammps.NeighOne = 8000
sim.export.lammps.UseTable2 = True
sim.export.lammps.InnerCutoff = 1.e-6
sim.export.lammps.NPairPotentialBins = 1000
sim.export.lammps.LammpsExec = 'lmp_omp'
sim.export.lammps.UseLangevin = True
sim.export.lammps.OMP_NumThread = __OMP_NumThread__
sim.export.lammps.TableInterpolationStyle = 'spline' # More robust than spline for highly CG-ed systems
sim.srel.optimizetrajlammps.LammpsDelTempFiles = False
sim.srel.optimizetrajlammps.UseLangevin = True
sim.export.omm.platformName = 'OpenCL' # or 'OpenCL' or 'GPU' or 'CUDA'
sim.export.omm.device = -1 #-1 is default, let openmm choose its own platform.
sim.export.omm.NPairPotentialKnots = 500 #number of points used to spline-interpolate the potential
sim.export.omm.InnerCutoff = 0.001 #0.001 is default. Note that a small value is not necessary, like in the lammps export, because the omm export interpolates to zero
sim.srel.optimizetrajomm.OpenMMStepsMin = 0 #number of steps to minimize structure, 0 is default
sim.export.omm.UseTabulated = True

"""---End of inputs---"""
        
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
                    LJGaussParams.update({(atom1,atom2): [BHOH_HOH, kappa, LJGDist0, Cut, LJGSigma, LJGEpsilon, 'LJGauss{}_{}'.format(atom1,atom2)]})
                    IsFixedLJGauss.update({(atom1,atom2): [True, True, FixedDist0, FixedSigma, FixedEpsilon, True]})
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
                PSplineParams.update({(atom1,atom2): [PSplineNKnot, Cut, NonbondEneSlopeInit, 'PSpline{}_{}'.format(atom1,atom2)]})
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

"""Create systems"""
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
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, Units = Units)

    Systems.append(Sys)
    NAtoms.append(Sys.NAtom)
    NMols.append(Sys.NMol)

""" Run MD """
if RunMD:
    import analysis
for i, Sys in enumerate(Systems):
    Sys.ForceField.SetParamString(open(ForceFieldFile, 'r').read())
    if StepScales:
        scaledStepsEquil = StepsEquil * StepScales[i]
        scaledStepsProd = StepsProd * StepScales[i]
    else:
        scaledStepsEquil = StepsEquil
        scaledStepsProd = StepsProd

    #make initial pdb
    top = sim.traj.pdb.PdbWrite('{}_init.pdb'.format(Sys.Name))
    top.AddAction(Sys.Int, StepFreq = 9)
#    Sys.Int.Run(10)

    NAtomsPerChain = len(MolTypesDicts[i]['PAA'])
    DOP = len(MolTypesDicts[i]['PAA'])
    NP = NMolsDicts[i]['PAA']
    DOP = len(MolTypesDicts[i]['PAA'])
    if UseSim and RunMD:
        fobj = open('{}_measures.dat'.format(Sys.Name), 'w')
        Sys.Measures.VerboseOutput(fobj = fobj, StepFreq=StepsStride)
        Trj = sim.traj.lammps.LammpsWrite("{}_traj.lammpstrj".format(Sys.Name))
        Trj.AddAction(Sys.Int, StepFreq = StepsStride)
        print "Now conducting warm-up...."
        Sys.Int.Run(scaledStepsEquil)
        #Sys.Measures.Reset()
        print "Now running production runs...."
        Sys.Int.Run(scaledStepsProd)
        print "timing:", Sys.Int.TimeElapsed
        print "\n"

    elif UseLammps:
        TrajFile = Sys.Name+'_traj.dcd'
        ThermoLog = Sys.Name+'_'+'lammps.log'
        top = Sys.Name+'.pdb'
        if RunMD:
            ret = sim.export.lammps.MakeLammpsTraj(Sys, DelTempFiles = False, Prefix = Sys.Name+'_', TrajFile = 'traj.dcd', WriteRestart = True, RestartFile = MDRestartFile,
                                                       Verbose = True, NStepsMin = StepsMin, NStepsEquil = scaledStepsEquil, NStepsProd = scaledStepsProd,
                                                       WriteFreq = StepsStride, CalcPress = False , OutputDCD = True, ReturnTraj = False, Nevery=500, Nrepeat=1, Nfreq=500,
                                                       ThermoOutput = "step pe temp press vol", ThermoStyle = "col")
            # Converts LAMMPS.data to .pdb with structure information
            u = mda.Universe(Sys.Name+'_'+'lammps.data')
            gr = u.atoms
            gr.write(top)
         
            """Analyze"""
            analysis.GetStats(TrajFile, top, NP, ThermoLog, DOP = DOP, NAtomsPerChain = NAtomsPerChain, StatsFName = 'AllStats.dat',
              RgDatName = 'RgTimeSeries', ReeDatName = 'ReeTimeSeries',RgStatOutName = 'RgReeStats', Ext='.dat',  
              fi = 'lammps', obs = ['PotEng', 'Temp', 'Press'], cols = None,
              res0Id = 0, stride = 1, autowarmup = True, warmup = 100)
            analysis.GetCompressibility(TrajFile, top, 298.15, unit = 'nonDim', lengthScale = 0.31, trajFmt = 'lmp')
    elif UseOMM:
        top = Sys.Name+'_'+'initial.pdb'
        TrajFile = Sys.Name+'_traj.dcd'
        ThermoLog = Sys.Name+'_'+'prodlog.txt'
        if RunMD:
            ret = sim.export.omm.MakeOpenMMTraj(Sys, DelTempFiles = False, Prefix = Sys.Name+'_', TrajFile = 'traj.dcd',
                                                       Verbose = True, NStepsMin = StepsMin, NStepsEquil = scaledStepsEquil, NStepsProd = scaledStepsProd,
                                                       WriteFreq = StepsStride)

            """Analyze"""
            analysis.GetStats(TrajFile, top, NP, ThermoLog, DOP = DOP, NAtomsPerChain = NAtomsPerChain, StatsFName = 'AllStats.dat',
              RgDatName = 'RgTimeSeries', ReeDatName = 'ReeTimeSeries',RgStatOutName = 'RgReeStats', Ext='.dat',
              fi = 'openmm', obs = [ 'Potential_Energy_(kJ/mole)',   'Kinetic_Energy_(kJ/mole)',   'Temperature_(K)',   'Box_Volume_(nm^3)',   'Density_(g/mL)'], cols = None,
              res0Id = 0, stride = 1, autowarmup = True, warmup = 100)
            analysis.GetCompressibility(TrajFile, top, 298.15, unit = 'nonDim', lengthScale = 0.31, trajFmt = 'omm')
    '''Do FEP'''
    if CalChemPot:

        if UseOMM:
            ThermoFile = ThermoLog
            PEColId = 1
            SetReRunRefState = SetReRunRefState
        elif UseLammps:
            ThermoFile = 'log_production'
            SetReRunRefState = True #always rerun lammps traj to make sure potential energy is in the correct unit
            PEColId = 1
        TopFile = top 

        os.mkdir(FEPDir)
        TrajFile =  os.path.join(cwd,TrajFile)
        TopFile = os.path.join(cwd,TopFile)
        ThermoFile = os.path.join(cwd,ThermoFile)
        ForceFieldFile = os.path.join(cwd,ForceFieldFile)
        os.chdir(FEPDir)
        
        dF, dF_stderr, muExcess, nFrames = fepcg.FEP(Sys, FEPMolNames, ThermoSlice, TrajSlice, Warmup, TrajFile, TopFile, ThermoFile, nDraw, nInsert, nDelete,
                              SetReRunRefState, PEColId,
                              BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams, ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, Units)

        s = 'Excess chemical potential calculation for {}'.format(FEPMolNames)
        s += '\nNumber of insertions: {}\nNumber of deletions: {}\nNumber of traj frames: {}'.format(nInsert,nDelete,nFrames)
        s += "\ndF: {0:4.4f} +/- {1:2.5f} kT".format(dF, dF_stderr)
        s += "\nBennett's Constant: {0:4.4f} kT".format(float(muExcess))
        print(s)
        f = open('chemicalPot.dat','w')
        f.write(s)
