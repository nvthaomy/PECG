#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:25:20 2019

@author: my
"""
import sim
import numpy as np
import forcefield
def CreateSystem(SysName, BoxL, UniqueCGatomTypes, MolNames, MolTypesDict, NMolsDict, charges, IsFixedCharge, Temp, Pres, IntParams, ForceFieldFile,
                              NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams, BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot, 
                              Units = sim.units.AtomicUnits,RgConstrain=False, MolIdRgs=[],IsFixedExtPot = {"UConst": True, "NPeriods":True}, StepsStride=1, RandomMul=0, RLength_dict={}, elements={}, nMonomers=0, L=[1.,1.,1.], a=0.31, polyL=None):

    print("\nCreate system {}".format(SysName))
    AtomTypes = {}
    MolTypes = []
    NMons = []
    IsCharged = False
    for AtomName in UniqueCGatomTypes:
        Charge = charges[AtomName]
        if Charge != 0:
            IsCharged = True
        if AtomName in elements.keys():
            Element = elements[AtomName]
        else:
            Element = None
        AtomType = sim.chem.AtomType(AtomName, Mass = 1., Charge = Charge, Element=Element)
        AtomTypes.update({AtomName:AtomType})
    print('AtomTypes in World: {}'.format(AtomTypes)) 
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
    print('MolTypes in World {}'.format(MolTypes)) 
    World = sim.chem.World(MolTypes, Dim = 3, Units = Units)

    #create bonds between monomer pairs
    for i,MolType in enumerate(MolTypes):
        NMon = NMons[i]        
        COMPos = [[0.,0.,0.]] # assume linear molecules
        for bond_index in range(0, NMon-1):
            if MolType.Name in RLength_dict.keys():
                print('Add rigid bond ', RLength_dict[MolType.Name][bond_index])
                MolType.Bond(bond_index, bond_index+1, RLength = RLength_dict[MolType.Name][bond_index])        
                COMPos.append(np.array(COMPos[-1]) + np.array([RLength_dict[MolType.Name][bond_index], 0., 0.]))
            else:
                MolType.Bond(bond_index, bond_index+1)
        if MolType.Name in RLength_dict.keys():
            COMPos = np.array(COMPos)
            MolType.COMPos = COMPos 
            print('COM Pos of {}: '.format(MolType.Name),COMPos)
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
    if not nMonomers:
        sim.system.positions.CubicLattice(Sys, Random = RandomMul)
    else:
        # center monomers at the center of an elongated box, assume elongated in the z direction
        import numpy as np
        Lz = L[2]
        Lx = L[0]
        Ly = L[1]
        if not polyL:
            polyRho = 11.5 #nm-3
            polyV = float(nMonomers)/polyRho/a**3 # minimum volume to pack polymers, convert to dimensionless units
            polyL = polyV/Lx/Ly * 1.5 # multiply by 1.5 to avoid overlap
        polyPos = np.random.rand(nMonomers,3) * np.array([Lx,Ly,polyL]) + np.array([0.,0.,Lz/2. - 0.5*polyL])
        Sys.Pos = np.random.rand(Sys.NAtom,3) * np.array([Lx,Ly,Lz])
        Sys.Pos[:nMonomers,:] = polyPos

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

