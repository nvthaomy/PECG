#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 13:03:23 2019

@author: my
"""
import sim
import os
import numpy as np

def CreateOptimizer(Sys, CGtraj, UseLammps, UseOMM, UseSim, StepsEquil, StepsProd, StepsStride, StepScale, UseWPenalty, ElecSys=None, RgConstrain=False,RgTars=[1.],measureRgs=None,recalc=False, LagMultList=np.zeros(500)):
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
    Opt = OptClass(Sys, Map, Beta = 1./Sys.TempSet, Traj = CGtraj, FilePrefix = '{}'.format(Sys.Name), ElecSys = ElecSys,
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
        Opt.AddPenalty("Virial", W, MeasureScale = 1./Sys.NAtom, Coef = 1.e-80, LagMult = LagMultList[0]) #HERE also need to scale the measure by 1/NAtom to be comparable to Srel

    # add Rg constraints for multiple species
    if RgConstrain:
        for i,RgTar in enumerate(RgTars):
            Opt.AddPenalty(measureRgs[i], RgTar, MeasureScale = 1., Coef = 1.e-80, LagMult = LagMultList[i+UseWPenalty])
    if recalc:
        Opt.CheckReady()
        Opt.StartIter = Opt.Iter
        SearchDirIter = 0
        Opt.Mode = "INIT"
        Opt.Backtracking = False
        Opt.UseMaxChange = False
        Opt.PrevLineSearch = False
        Opt.CGRestartCount = 0
        Opt.Initdx()

        Opt.InitConstraints(Verbose = Opt.Verbose)
        Opt.SetParam(Opt.Param)
        if Opt.Verbose:
            print Opt.Output0()
        Opt.OutputTarHistFile()
        Opt.ReweightTar = False #i.e. make a new model trajectory
    return Opt

def RunOpt(Opts, Weights, Prefix, UseWPenalty, MaxIter, SteepestIter, RgConstrain, StageCoefs=[1e8, 1e10, 1e12],UpdateMode = 0, NMaxStage = 100):

    Optimizer = sim.srel.OptimizeMultiTrajClass(Opts, Weights=Weights)
    Optimizer.FilePrefix = ("{}".format(Prefix))
    if not UseWPenalty and not RgConstrain:
        Optimizer.RunConjugateGradient(MaxIter=MaxIter, SteepestIter=SteepestIter)
    else:
        Optimizer.RunStages(StageCoefs = StageCoefs, UpdateMode = UpdateMode) #, NMaxStage = NMaxStage)
def recalc(Opts,Prefix):
    import time
    StartTime = time.time()
    #Opt = Opts[0] 

    Opt = sim.srel.OptimizeMultiTrajClass(Opts)
    Opt.FilePrefix = ("{}".format(Prefix))

    print('===== Begin recalculating trajectory and parametric derivatives =====')
    if not Opt.ReweightTar:
        Opt.UpdateModTraj()
        Opt.OutputModHistFile()
        Opt.OutputPlot()

    NotFixed = np.logical_not(Opt.Fixed)
    Opt.CalcObj(CalcDeriv = True)

    print("=== Hessian/DDSrel ===")
    np.savetxt("{}_Hessian.dat".format(Prefix),Opt.DDObj)

    Ignore = Opt.Fixed.copy()
    for i in Opt.Constrained:
        Ignore[i] = True

    for (i, ThisDObj) in enumerate(Opt.DDObj):
        if all(Opt.DDObj[i,:] == 0):
            Ignore[i] = True

    Use = np.logical_not(Ignore)

    DObj = Opt.DObj[Use]
    DDObj = Opt.DDObj[np.ix_(Use, Use)] + Opt.HessianPad * np.identity(len(DObj))

    print('Used variables: {}'.format(Use))
    print('Non-zero Hessian:\n{}'.format(DDObj))
    np.savetxt('{}_Hessian_masked.dat'.format(Prefix),DDObj)

    print("=== Gradient/DSrel ===")
    np.savetxt('{}_grad.dat'.format(Prefix),Opt.DObj)
    print('Non-zero Grad:\n{}'.format(DObj))
    np.savetxt('{}_grad_masked.dat'.format(Prefix),DObj)

    print("=== Misc. ===")
    print("Total recalculation runtime: {}min".format( (time.time()-StartTime)/60. ))

