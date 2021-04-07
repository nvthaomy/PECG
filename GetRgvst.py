#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:18:47 2021

@author: nvthaomy
"""
import numpy as np
import sys

file = open(sys.argv[1],'r')

outfile = 'Rg_vs_t.txt'
nMeasure = 1

'============================='
iter = []
deltaRg = {}
lines = file.readlines()
readRg = False
for i in range(nMeasure):
    deltaRg.update({'m{}'.format(i):[]})
for j,l in enumerate(lines):
    if 'ITERATION' in l:
        iter.append(float(l.split()[1]))
        readRg = True
        i = 0
    if 'Delta' in l:
        deltaRg['m{}'.format(i)].append(float(l.split()[-1]))
        i+= 1
n = min(len(iter),len(deltaRg['m{}'.format(nMeasure-1)]))
deltaRg_tmp = deltaRg
deltaRg = {}
for k,v in deltaRg_tmp.items():
    deltaRg.update({k:v[:n]})
keys=list(deltaRg.keys())
deltaRg=list(deltaRg.values())

dat = np.transpose(np.vstack((iter[:n],*deltaRg)))
np.savetxt(outfile,dat,delimiter=',',header='iter {}'.format(' '.join(keys)),fmt=['%i','%.5e'])
