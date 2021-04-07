#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:18:47 2021

@author: nvthaomy
"""
import numpy as np
import sys

file = open(sys.argv[1],'r')

outfile = 'SimDerivatives.txt'
Pname = 'BondA_A:Dist0'
header = '=== R0 '

'============================='
coef = []
lagmul = []
P = []
Rg = []
DSrelDl = []
DObjDl = []
DBiasDl = []
DAvgDl = []

lines = file.readlines()
readP = False
for j,l in enumerate(lines):
    if 'COEF =' in l:
        coef.append(float(l.split()[-1]))
    if 'LAGMULT' in l:
        lagmul.append(float(l.split()[-1]))
    if 'New Calculate Avg value' in l:
        Rg.append(float(l.split()[-1]))
    if 'INITIAL PARAMETERS' in l:
        readP = True
    if readP and Pname in l:
        P.append(float(l.split()[-1]))
        readP = False
    if 'DObj' in l and '{}'.format(Pname) in l:
        DObjDl.append(float(l.split()[-1]))
        DSrelDl.append(float(l.split()[-2]))
        DBiasDl.append(DObjDl[-1]-DSrelDl[-1])
    if 'DAvg' in l and '{}'.format(Pname) in l:
        DAvgDl.append(float(l.split()[-1]))

n = len(DAvgDl)
dat = np.transpose(np.vstack((P[:n],DSrelDl[:n],DObjDl[:n],DBiasDl[:n],DAvgDl[:n],Rg[:n],coef[:n],lagmul[:n])))     
np.savetxt(outfile,dat,delimiter=',',header='{} DSrelDl DObjDl DBiasDl DRgDl Rg Coef LagMult'.format(Pname),fmt='%.8e')
