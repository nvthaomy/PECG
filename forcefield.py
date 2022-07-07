#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:17:26 2019

@author: my
"""                        
import sim
print(sim)

def CreateForceField(Sys, IsCharged, AtomTypes, NGaussDicts, LJGaussParams, IsFixedLJGauss, SmearedCoulParams, EwaldParams,
                              BondParams, IsFixedBond, PSplineParams, UseLJGauss, ExtPot,IsFixedExtPot):
    """BondParams: (atom1,atom2):[Dist0,FConst,Label]
       IsFixedBond: (atom1,atom2):[Dist0,FConst,Label], type = boolean
       LJGaussParams: (atom1,atom2):[B, kappa, Dist0, Cut, Sigma, Epsilon, Label ]
       IsFixedLJGauss: (atom1,atom2):[B, kappa, Dist0, Cut, Sigma, Epsilon, Label ] type = boolean
       SmearedCoulParams: (atom1,atom2): [BornA, Cut, Shift, FixedCoef, FixedBornA, Label]
       EwaldParams = {'ExcludeBondOrd': 3, 'Cut': Cut, 'Shift': True, 'Label': 'EW'}
       ExtPot = {"UConst": 0, "NPeriods": 1, "PlaneAxis": 2, "PlaneLoc": 0., "AtomTypes":[]}
       IsFixedExtPot = {"UConst": boolean, "NPeriods": boolean}"""
   
    print("\nCreating forcefield for {}".format(Sys.Name))
    ForceField = []
    try:
        if ExtPot["UConst"] > 0:
            UseExternal = True
        else:
            UseExternal = False
    except: 
        for P in ExtPot:
            if P["UConst"] > 0:
                UseExternal = True
                break
            else: 
                UseExternal = False
    #Bond    
    for (atom1name,atom2name), params in sorted(BondParams.items()):
        print('Adding {}'.format(params[2]))
        if not atom1name == 'All' and not atom2name == 'All':
            atom1 = AtomTypes[atom1name]
            atom2 = AtomTypes[atom2name]
            Filter = sim.atomselect.PolyFilter(Filters = [atom1, atom2], Bonded = True)
        else:
            Filter = sim.atomselect.BondPairs
        fixed = IsFixedBond[(atom1name,atom2name)]  
        P = sim.potential.Bond(Sys, Filter = Filter, Dist0 = params[0], FConst = params[1], Label = params[2])
        #fixing params
        P.Param.Dist0.Fixed = fixed[0]
        P.Param.FConst.Fixed = fixed[1]
        ForceField.append(P)
    if not UseLJGauss:
        #Spline Pair
        for (atom1name,atom2name), params in sorted(PSplineParams.items()):
            print('Adding {}\nCutoff {}'.format(params[3],params[1]))
            atom1 = AtomTypes[atom1name]
            atom2 = AtomTypes[atom2name]
            Filter = sim.atomselect.PolyFilter(Filters = [atom1, atom2])
            P = sim.potential.PairSpline(Sys, Filter = Filter, Cut = params[1], Fixed = params[4],
                                           NKnot = params[0], Label = params[3],
                                           NonbondEneSlopeInit = params[2])
            ForceField.append(P)
    elif UseLJGauss:
        #Gaussian Pair
        for (atom1name,atom2name), params in sorted(LJGaussParams.items()):
            print('Adding {}\nCutoff {}'.format(params[-1],params[3]))
            if 'MOL_' in atom1name:
                MolType = [x for x in Sys.World if x.Name == ''.join(atom1name.split('MOL_'))][0]
                atom1 = sim.atomselect.Filter(MolTypes = MolType)
            elif 'ELE_' in atom1name:
                atom1 = sim.atomselect.Filter(Element = ''.join(atom1name.split('ELE_')))
            else:
                try:
                    atom1 = AtomTypes[atom1name]
                except:
                    atom1 = [AtomTypes[x] for x in atom1name]
            
            if 'MOL_' in atom2name:
                MolType = [x for x in Sys.World if x.Name == ''.join(atom2name.split('MOL_'))][0]
                atom2 = sim.atomselect.Filter(MolTypes = MolType)
            elif 'ELE_' in atom2name:
                atom2 = sim.atomselect.Filter(Element = ''.join(atom2name.split('ELE_')))
            else:
                try:
                    atom2 = AtomTypes[atom2name]
                except:
                    atom2 = [AtomTypes[x] for x in atom2name]

            fixed = IsFixedLJGauss[(atom1name,atom2name)]
            print('atom 1: {}'.format(atom1))
            print('atom 2: {}'.format(atom2))
            Filter = sim.atomselect.PolyFilter(Filters = [atom1, atom2])
            print('Filter ', Filter)
            if ('All','All') in NGaussDicts.keys():
                NGauss = NGaussDicts[('All','All')]
            else:
                try:
                    NGauss = NGaussDicts[(atom1name,atom2name)]
                except:
                    NGauss = NGaussDicts[(atom2name,atom1name)]
            #adding multiple Gaussians
            for i in range(NGauss): 
                Label = params[-1] + '{}'.format(i)
                P = sim.potential.LJGaussian(Sys, Filter = Filter, B = params[0], Kappa = params[1],
                             Dist0 = params[2], Cut = params[3], Sigma = params[4], Epsilon = params[5], 
                             Label = Label)
                P.Param.B.Fixed = fixed[0]
                #set bound of B, repulsive if i is even, attractive if i is odd
                if i%2 == 0:
                    P.Param.B.Min = 0.
                    #if  params[0] < 0.:
                    #    P.Param.B = -params[0]
                else:
                    P.Param.B.Max = 0.
                    if params[0] > 0.:
                        P.Param.B = -params[0]
                P.Param.Kappa.Fixed = fixed[1]
                P.Param.Dist0.Fixed = fixed[2]
                P.Param.Sigma.Fixed = fixed[4]
                P.Param.Epsilon.Fixed = fixed[5]
                ForceField.append(P)

    if IsCharged:
        #Ewald
        print('Adding Ewald')
        P = sim.potential.Ewald(Sys, ExcludeBondOrd = EwaldParams['ExcludeBondOrd'] , 
                         Cut = EwaldParams['Cut'], Shift = EwaldParams['Shift'], Label = EwaldParams['Label'], Coef = EwaldParams['Coef'],# EwaldNAtom = EwaldParams['EwaldNAtom'],
                         FixedCoef = EwaldParams['FixedCoef']) #, KMax = EwaldParams['KMax'])
#        print('Number of atoms to loop through in Ewald {}'.format(P.EwaldNAtom))
        ForceField.append(P)
    
        #Smeared Coulomb
        for (atom1name,atom2name), params in sorted(SmearedCoulParams.items()):
            print('Adding {}'.format(params[5]))
            atom1 = AtomTypes[atom1name]
            atom2 = AtomTypes[atom2name]
            Filter = sim.atomselect.PolyFilter(Filters = [atom1, atom2])
            P = sim.potential.SmearedCoulombEwCorr(Sys, Filter = Filter, BornA = params[0], Cut = params[1], Shift = params[2],
                                                FixedCoef = params[3], FixedBornA = params[4], Label = params[5], Coef = params[6])
            ForceField.append(P)
        
    if UseExternal:
        """applies external field to just species included in FilterExt"""
        if isinstance(ExtPot,list):
            for P in ExtPot: 
                AtomTypesInExt = AtomTypes[P['AtomTypes']]
                FilterExt = sim.atomselect.PolyFilter(Filters = [AtomTypesInExt])
                print("Using external sinusoid with UConst {} on {}".format(P["UConst"],P['AtomTypes'])) 
                P0 = sim.potential.ExternalSinusoid(Sys, Filter=FilterExt, UConst=P["UConst"], NPeriods=P["NPeriods"], 
                                           PlaneAxis=P["PlaneAxis"], PlaneLoc=P["PlaneLoc"], Label=P["Label"])
                P0.UConst.Fixed = IsFixedExtPot['UConst']  
                P0.NPeriods.Fixed = IsFixedExtPot['NPeriods']
                ForceField.append(P0)
        else:
            print("Using external sinusoid with UConst {} on {}".format(ExtPot["UConst"],ExtPot['AtomTypes']))    
            AtomTypesInExt = AtomTypes[ExtPot['AtomTypes']]
            FilterExt = sim.atomselect.PolyFilter(Filters =[AtomTypesInExt]) 
            P0 = sim.potential.ExternalSinusoid(Sys, Filter=FilterExt, UConst=ExtPot["UConst"], NPeriods=ExtPot["NPeriods"],     
                          PlaneAxis=ExtPot["PlaneAxis"], PlaneLoc=ExtPot["PlaneLoc"], Label=ExtPot["Label"])

            P0.UConst.Fixed = IsFixedExtPot['UConst']
            P0.NPeriods.Fixed = IsFixedExtPot['NPeriods']
            ForceField.append(P0)
    return ForceField
    
    
    
