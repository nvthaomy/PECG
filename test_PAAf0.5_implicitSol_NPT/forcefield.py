#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:17:26 2019

@author: my
"""                        
import sim
print(sim)

def CreateForceField(Sys, IsCharged, AtomTypes, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams,
                              BondParams, IsFixedBond, ExtPot):
    """BondParams: (atom1,atom2):[Dist0,FConst,Label]
       IsFixedBond: (atom1,atom2):[Dist0,FConst,Label], type = boolean
       LJGaussParams: (atom1,atom2):[B, kappa, Dist0, Cut, Sigma, Epsilon, Label ]
       IsFixedLJGauss: (atom1,atom2):[B, kappa, Dist0, Cut, Sigma, Epsilon, Label ] type = boolean
       SmearedCoulParams: (atom1,atom2): [BornA, Cut, Shift, FixedCoef, FixedBornA, Label]
       EwaldParams = {'ExcludeBondOrd': 3, 'Cut': Cut, 'Shift': True, 'Label': 'EW'}
       ExtPot = {"UConst": 0, "NPeriods": 1, "PlaneAxis": 2, "PlaneLoc": 0., "AtomTypes":[]}"""
   
    print("\nCreating forcefield for {}".format(Sys.Name))
    ForceField = []
    if ExtPot["UConst"] > 0:
        UseExternal = True
    else:
        UseExternal = False
    #Bond    
    for (atom1name,atom2name), params in sorted(BondParams.items()):
        print('Adding {}'.format(params[2]))
        atom1 = AtomTypes[atom1name]
        atom2 = AtomTypes[atom2name]
        fixed = IsFixedBond[(atom1name,atom2name)]
        Filter = sim.atomselect.PolyFilter(Filters = [atom1, atom2], Bonded = True)
        P = sim.potential.Bond(Sys, Filter = Filter, Dist0 = params[0], FConst = params[1], Label = params[2])
        #fixing params
        P.Param.Dist0.Fixed = fixed[0]
        P.Param.FConst.Fixed = fixed[1]
        ForceField.append(P)
    #Gaussian Pair
    for (atom1name,atom2name), params in sorted(LJGaussParams.items()):
        print('Adding {}\nCutoff {}'.format(params[6],params[3]))
        atom1 = AtomTypes[atom1name]
        atom2 = AtomTypes[atom2name]
        fixed = IsFixedLJGauss[(atom1name,atom2name)]
        Filter = sim.atomselect.PolyFilter(Filters = [atom1, atom2])
        P = sim.potential.LJGaussian(Sys, Filter = Filter, B = params[0], Kappa = params[1],
                             Dist0 = params[2], Cut = params[3], Sigma = params[4], Epsilon = params[5], 
                             Label = params[6])
        P.Param.B.Fixed = fixed[0]
        P.Param.Kappa.Fixed = fixed[1]
        P.Param.Dist0.Fixed = fixed[2]
        P.Param.Sigma.Fixed = fixed[4]
        P.Param.Epsilon.Fixed = fixed[5]
        ForceField.append(P)

    if IsCharged:
        #Ewald
        print('Adding Ewald')
        P = sim.potential.Ewald(Sys, ExcludeBondOrd = EwaldParams['ExcludeBondOrd'] , 
                         Cut = EwaldParams['Cut'], Shift = EwaldParams['Shift'], Label = EwaldParams['Label'])
        ForceField.append(P)
    
        #Smeared Coulomb
        for (atom1name,atom2name), params in sorted(SmearedCoulParams.items()):
            print('Adding {}'.format(params[5]))
            atom1 = AtomTypes[atom1name]
            atom2 = AtomTypes[atom2name]
            Filter = sim.atomselect.PolyFilter(Filters = [atom1, atom2])
            P = sim.potential.SmearedCoulombEwCorr(Sys, Filter = Filter, BornA = params[0], Cut = params[1], Shift = params[2],
                                                FixedCoef = params[3], FixedBornA = params[4], Label = params[5])
            ForceField.append(P)
        
    if UseExternal:
        """applies external field to just species included in FilterExt"""
        print("Using external sinusoid with UConst {}".format(ExtPot["UConst"]))
        AtomTypesInExt = []
        for AtomName in ExtPot['AtomTypes']:
            AtomTypesInExt.append(AtomTypes[AtomName])
        FilterExt = sim.atomselect.PolyFilter(AtomTypesInExt)
        P = sim.potential.ExternalSinusoid(Sys, Filter=FilterExt, UConst=ExtPot["UConst"], NPeriods=ExtPot["NPeriods"], 
                                           PlaneAxis=ExtPot["PlaneAxis"], PlaneLoc=ExtPot["PlaneLoc"], Label="ExtSin")
        ForceField.append(P)
    return ForceField
    
    
    
