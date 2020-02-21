#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:25:20 2019
Set up systems for chemical potential calculation from insertion and deletion
Tested for traj outputted by openmm

TO DO:
	need to shift trajectory by lammps?

@author: my
"""
import sim
import forcefield

def CreateSystem1(SysName, FEPMolNames, BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams, ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, 
                              Units):
    """Create system of additional residue for FEP
    Adding this residue to the end of World object"""

    print("\nCreate system {}".format(SysName))
    AtomTypes = {}
    MolTypes = []
    NMons = []
    IsCharged = False
    for AtomName in UniqueCGatomTypes:
        Charge = charges[AtomName]
        if Charge != 0:
            IsCharged = True
        AtomType = sim.chem.AtomType(AtomName, Mass = 1., Charge = Charge)
        AtomTypes.update({AtomName:AtomType})
    
    #need to create MolType in a right sequence to the trajectory
    for MolName in MolNames:  
        if NMolsDict[MolName] > 0: 
            AtomsInMol = []
            AtomNames = MolTypesDict[MolName]
            NMons.append(len(AtomNames))
            for AtomName in AtomNames:
                AtomsInMol.append(AtomTypes[AtomName])
            MolType = sim.chem.MolType(MolName, AtomsInMol)
            MolTypes.append(MolType)
    World = sim.chem.World(MolTypes, Dim = 3, Units = Units)

    #create bonds between monomer pairs
    for i,MolType in enumerate(MolTypes):
        NMon = NMons[i]        
        for bond_index in range(0, NMon-1):
            MolType.Bond(bond_index, bond_index+1)

    # make system    
    Sys = sim.system.System(World, Name = SysName)
    Sys.BoxL = BoxL
    for i, MolType in enumerate(MolTypes): 
        NMol = NMolsDict[MolType.Name]
        print("Adding {} {} molecules to system".format(NMol, MolType.Name))
        for j in range(NMol):
            Sys += MolType.New()
    print("Adding one molecule for molecule names: {}".format(FEPMolNames))
    for MolName in FEPMolNames:
        for MolType in MolTypes:
            if MolType.Name == MolName:
                Sys += MolType.New()
    print('{} molecules in system1'.format(Sys.NMol))
    Sys.ForceField.Globals.Charge.Fixed = IsFixedCharge
    
    # add forcefield 
    
    ForceField = forcefield.CreateForceField(Sys, IsCharged, AtomTypes, NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams,
                              BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot)
                                
    Sys.ForceField.extend(ForceField)

    ''' Now setting initial system optimizations. '''
    if ForceFieldFile: 
        with open(ForceFieldFile, 'r') as of: s = of.read()
        Sys.ForceField.SetParamString(s)      
        
    #set up the histograms
    for P in Sys.ForceField:
        P.Arg.SetupHist(NBin = 10000, ReportNBin=100)
    # lock and load
    Sys.Load()


    #initial positions and velocities
    sim.system.positions.CubicLattice(Sys)
    sim.system.velocities.Canonical(Sys, Temp = Temp)
    Sys.TempSet = Temp
    Sys.PresSet = Pres
    
    #configure integrator
    Int = Sys.Int
    Int.Method = Int.Methods.VVIntegrate        
    Int.Method.TimeStep = IntParams['TimeStep']

    if IntParams['ensemble'] == 'NVT':
        Int.Method.Thermostat = Int.Method.ThermostatLangevin
        Int.Method.LangevinGamma = IntParams['LangevinGamma']
    elif IntParams['ensemble'] == 'NPT':
        Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
        Int.Method.Barostat = Int.Method.BarostatMonteCarlo  
    
    return Sys

def CreateSystem2(SysName, FEPMolNames, BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams, ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot,
                              Units):
    """Create system of one less residue for FEP"""

    print("\nCreate system {}".format(SysName))
    AtomTypes = {}
    MolTypes = []
    NMons = []
    IsCharged = False
    for AtomName in UniqueCGatomTypes:
        Charge = charges[AtomName]
        if Charge != 0:
            IsCharged = True
        AtomType = sim.chem.AtomType(AtomName, Mass = 1., Charge = Charge)
        AtomTypes.update({AtomName:AtomType})

    #need to create MolType in a right sequence to the trajectory
    for MolName in MolNames:
        if NMolsDict[MolName] > 0:
            AtomsInMol = []
            AtomNames = MolTypesDict[MolName]
            NMons.append(len(AtomNames))
            for AtomName in AtomNames:
                AtomsInMol.append(AtomTypes[AtomName])
            MolType = sim.chem.MolType(MolName, AtomsInMol)
            MolTypes.append(MolType)
    World = sim.chem.World(MolTypes, Dim = 3, Units = Units)

    #create bonds between monomer pairs
    for i,MolType in enumerate(MolTypes):
        NMon = NMons[i]
        for bond_index in range(0, NMon-1):
            MolType.Bond(bond_index, bond_index+1)

    # make system    
    Sys = sim.system.System(World, Name = SysName)
    Sys.BoxL = BoxL
    for i, MolType in enumerate(MolTypes):
        if not MolType.Name in FEPMolNames:
            NMol = NMolsDict[MolType.Name]
        else:
            NMol = NMolsDict[MolType.Name] - 1
        print("Adding {} {} molecules to system".format(NMol, MolType.Name))
        for j in range(NMol):
            Sys += MolType.New()
    print('{} molecules in system2'.format(Sys.NMol))

    Sys.ForceField.Globals.Charge.Fixed = IsFixedCharge

    # add forcefield 

    ForceField = forcefield.CreateForceField(Sys, IsCharged, AtomTypes, NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams,
                              BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot)

    Sys.ForceField.extend(ForceField)

    ''' Now setting initial system optimizations. '''
    if ForceFieldFile:
        with open(ForceFieldFile, 'r') as of: s = of.read()
        Sys.ForceField.SetParamString(s)

    #set up the histograms
    for P in Sys.ForceField:
        P.Arg.SetupHist(NBin = 10000, ReportNBin=100)
    # lock and load
    Sys.Load()


    #initial positions and velocities
    sim.system.positions.CubicLattice(Sys)
    sim.system.velocities.Canonical(Sys, Temp = Temp)
    Sys.TempSet = Temp
    Sys.PresSet = Pres

    #configure integrator
    Int = Sys.Int
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.TimeStep = IntParams['TimeStep']

    if IntParams['ensemble'] == 'NVT':
        Int.Method.Thermostat = Int.Method.ThermostatLangevin
        Int.Method.LangevinGamma = IntParams['LangevinGamma']
    elif IntParams['ensemble'] == 'NPT':
        Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
        Int.Method.Barostat = Int.Method.BarostatMonteCarlo

    return Sys

def CreateFEPSystems(SysName, FEPMolNames,*args): 
    """create Sys1 (N+1 molecules) and Sys2 (N-1 molecules)"""
    Sys1 = CreateSystem1(SysName+'_1', FEPMolNames, *args)
    Sys2 = CreateSystem2(SysName+'_2', FEPMolNames, *args)
    return Sys1, Sys2

def FEP(Sys, FEPMolNames, ThermoSlice, TrajSlice, Warmup, TrajFile, TopFile, ThermoFile, nDraw, nInsert, nDelete, 
        SetReRunRefState, PEColId, *args):
                              # *args: BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams, ForceFieldFile,
                              #NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, Units
                              
    import mdtraj as md
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
    import FEP_Module_CG as FEP
    import numpy as np 
 
    Sys1, Sys2 = CreateFEPSystems(Sys.Name, FEPMolNames, *args)
    #Create OMM objects
    Verbose=False
    system0,top0,mdtrajtop0 = sim.export.omm.CreateOMMSys(Sys,Verbose)
    system1,top1,mdtrajtop1 = sim.export.omm.CreateOMMSys(Sys1,Verbose)
    system2,top2,mdtrajtop2 = sim.export.omm.CreateOMMSys(Sys2,Verbose)
    simulation0,simOptions0 = sim.export.omm.CreateOMMSimulation(Sys,system0,top0,Sys.Name+'_',Verbose=Verbose)
    simulation1,simOptions1 = sim.export.omm.CreateOMMSimulation(Sys1,system1,top1,Sys1.Name+'_',Verbose=Verbose)
    simulation2,simOptions2 = sim.export.omm.CreateOMMSimulation(Sys2,system2,top2,Sys2.Name+'_',Verbose=Verbose)
    temp = simOptions0['temp'].value_in_unit(kelvin)

    traj = md.load(TrajFile, top=TopFile, stride=TrajSlice)
    traj = traj[Warmup::]
    PotEneData = np.loadtxt(ThermoFile,delimiter='\t')[::ThermoSlice, PEColId]
    PotEneData = PotEneData[Warmup::]

    print('Number of frames in trajectory: {}'.format(len(traj)))

    if len(traj) != len(PotEneData):
        Exception("WARNING: The number of entries in the trajectory and thermo. file do not match!")

    print('\n=== Calculate the Chemical Potential (Excess) ===')
    # Create an FEP object
    FEP_Object = FEP.FEP('OB','particle',[traj],[simulation0,simulation1,simulation2],[PotEneData])

    # Set a variaty of object attributes
    FEP_Object.SetMethod('OB') # optimal-Bennett's method
    FEP_Object.SetDerivative('particle') # particle derivative
    FEP_Object.SetResidueName(FEPMolNames) # the residue type to be inserted/deleted;
    FEP_Object.SetNumberOfDrawsPerFrame(nDraw) # How many molecule samples to draw from on each frame to build insertion library.
    FEP_Object.SetNumberInsertionsPerFrame(nInsert)
    FEP_Object.SetNumberDeletionsPerFrame(nDelete)
    FEP_Object.SetReRunRefState(SetReRunRefState)
    FEP_Object.SetTemperature(temp) 

    # Names and Numbers of each residue, and the number of atoms in the residue for each trajectory
    print('Names and Number of each residue in each trajectory:')
    print(FEP_Object.GetResidueNames())

    # Perform the reweighting
    FEP_Object.Reweight()

    # Calculate the free energy difference
    FEP_Object.SetBennettsConstant(-33.8) # Initial guess, KJ/mole
    FEP_Object.Optimal_Bennetts() # Run optimal-Bennetts

    return FEP_Object.dF,FEP_Object.dF_stderr, FEP_Object.BennettsConstant, len(traj)

