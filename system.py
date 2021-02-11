#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:25:20 2019

@author: my
"""
import sim
import forcefield
def CreateSystem(SysName, BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams, ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, 
                              Units = sim.units.AtomicUnits,RgConstrain=False, MolIdRgs=[],IsFixedExtPot = {"UConst": True, "NPeriods":True}, StepsStride=1):

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
    
    #make dictionary to keep track of index use to set NMol
    nextNMolId = {}
    for MolName in MolNames:
        if not MolName in nextNMolId.keys():
            nextNMolId.update({MolName: 0})
    #need to create MolType in a right sequence to the trajectory
    for MolName in MolNames:  
        if (NMolsDict[MolName] > 0 and not isinstance(NMolsDict[MolName],list)) or (isinstance(NMolsDict[MolName],list)): 
            if isinstance(NMolsDict[MolName],list) and (not any(NMolsDict[MolName])or len(NMolsDict[MolName])==0):
                Exception('Does not handle value of 0 if NMol is list. If is empty list, set to 0')
            AtomsInMol = []
            AtomNames = MolTypesDict[MolName]
            NMons.append(len(AtomNames))
            for AtomName in AtomNames:
                AtomsInMol.append(AtomTypes[AtomName])
            MolType = sim.chem.MolType(MolName, AtomsInMol)
            MolTypes.append(MolType)
    print('MolTypes {}'.format(MolTypes)) 
    World = sim.chem.World(MolTypes, Dim = 3, Units = Units)

    #create bonds between monomer pairs
    for i,MolType in enumerate(MolTypes):
        NMon = NMons[i]        
        for bond_index in range(0, NMon-1):
            MolType.Bond(bond_index, bond_index+1)

    # make system and add molecules in same sequence as MolNames 
    Sys = sim.system.System(World, Name = SysName)
    Sys.BoxL = BoxL
    for i, MolType in enumerate(MolTypes): 
        NMol = NMolsDict[MolType.Name]
        if isinstance(NMol,list):
            NMol = NMol[nextNMolId[MolType.Name]]
            nextNMolId[MolType.Name] += 1
        print("Adding {} {} molecules to system".format(NMol, MolType.Name))
        for j in range(NMol):
            Sys += MolType.New()
    Sys.ForceField.Globals.Charge.Fixed = IsFixedCharge
    
    # add forcefield 
    
    ForceField = forcefield.CreateForceField(Sys, IsCharged, AtomTypes, NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams,
                              BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, IsFixedExtPot)
                                
    Sys.ForceField.extend(ForceField)

    ''' Now setting initial system optimizations. '''
    if ForceFieldFile: 
        with open(ForceFieldFile, 'r') as of: s = of.read()
        Sys.ForceField.SetParamString(s)      
        
    if RgConstrain == True:
        measureRgs = []
        try:
            for i,MolIdRg in enumerate(MolIdRgs):
                print('Adding RgEnsemble (new) measurement for molecules of indices {}'.format(MolIdRg))
                measureRg = sim.measure.rg.Rg(Sys, StepFreq = StepsStride, MolIndices=MolIdRg)
                Sys.Measures.append(measureRg)
                measureRgs.append(measureRg)
        except:
            print('Failed adding Rg (new) measurement')
            try:    
                for i,MolIdRg in enumerate(MolIdRgs):
                    print('Adding RgEnsemble (old) measurement for molecules of indices {}'.format(MolIdRg))
                    measureRg = sim.measure.rg.RgEnsemble(Sys, StepFreq = StepsStride, MolIndices=MolIdRg) 
                    Sys.Measures.append(measureRg)
                    measureRgs.append(measureRg)
            except:
                print('Failed adding RgEnsemble measurement')
    else:
        measureRgs = []
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
    
    return Sys,measureRgs 

