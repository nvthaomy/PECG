#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:18:47 2021

@author: nvthaomy
"""
import numpy as np
import sys

file = open(sys.argv[1],'r')

outfile = 'SimFobj_vs_k_LagOn.txt'
Pname = 'BondA_A:Dist0'

'============================='
coef = []
lagmul = []
P = []
Rg = []
lines = file.readlines()
readP = False
for j,l in enumerate(lines):
    if 'COEF =' in l:
        coef.append(float(l.split()[-1]))
    if 'LAGMULT' in l:
        lagmul.append(float(l.split()[-1]))
    if 'New Calculate Avg value' in l:
        ind = j
    if 'FINAL VALUES' in l:
        readP = True
        Rg.append(float(lines[ind].split()[-1]))
    if readP and Pname in l:
        P.append(float(l.split()[-1]))
        readP = False
n = min(len(P),len(coef),len(lagmul),len(Rg))
dat = np.transpose(np.vstack((coef[:n],lagmul[:n],P[:n],Rg[:n])))     
np.savetxt(outfile,dat,delimiter=',',header='coef lag_mult {} Rg'.format(Pname))
