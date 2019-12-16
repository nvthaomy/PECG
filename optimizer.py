#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 13:03:23 2019

@author: my
"""
import sim
import os
import numpy as np

def CreateOptimizer(Sys, CGtraj, UseLammps, UseOMM, UseSim, StepsEquil, StepsProd, StepsStride, StepScale, UseWPenalty):
    # Perform atom mapping for specific system
    Map = sim.atommap.PosMap()
    print(Sys.Name)
    print('NMol: {}'.format(Sys.NMol))
    print('NAtom: {}'.format(Sys.NAtom))
    print('NDOF: {}'.format(Sys.NDOF))
    for (i, a) in enumerate(Sys.Atom):
        Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]   

    ''' Setup Optimizers '''
    if UseLammps:
        OptClass = sim.srel.optimizetrajlammps.OptimizeTrajLammpsClass
    elif UseOMM:
        OptClass = sim.srel.optimizetrajomm.OptimizeTrajOpenMMClass
    elif UseSim:
        OptClass = sim.srel.optimizetraj.OptimizeTrajClass
        fobj = open('tmp_measures{}.dat'.format(Sys.Name), 'w')
        Sys.Measures.VerboseOutput(fobj = fobj, StepFreq = StepsStride)
        Int = Sys.Int
        Trj = sim.traj.lammps.LammpsWrite("tmp_traj_{}.lammpstrj".format(Sys.Name))
        Trj.AddAction(Int, StepFreq = StepsStride)
              
    Sys.ScaleBox(Sys.BoxL) # scale the system by the box
    Opt = OptClass(Sys, Map, Beta = 1./Sys.TempSet, Traj = CGtraj, FilePrefix = '{}'.format(Sys.Name),
                        SaveLoadArgData = True, TempFileDir = os.getcwd(), UseTarHists=False)                        
    Opt.ConstrainNeutralCharge()
    # Set run times for optimization objects.
    # Useful for dilute systems versus concentration systems.    
    if StepScale:
        temp_StepsEquil = StepsEquil * StepScale
        temp_StepsProd  = StepsProd * StepScale
    else:
        temp_StepsEquil  = StepsEquil
        temp_StepsProd   = StepsProd
    
    Opt.StepsEquil  = temp_StepsEquil
    Opt.StepsProd   = temp_StepsProd
    Opt.StepsStride = StepsStride    
    # NEED TO CHECK HOW THE VIRIAL IS COMPUTED BELOW, WAS ONLY USED FOR GAUSSIAN FLUID
    if UseWPenalty == True:
        Volume = np.prod(Sys.BoxL)
        W = Sys.NDOF - 3*Sys.Pres*Volume
        Opt.AddPenalty("Virial", W, MeasureScale = 1./Sys.NAtom, Coef = 1.e-80) #HERE also need to scale the measure by 1/NAtom to be comparable to Srel
    return Opt
