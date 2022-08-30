# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:29:13 2019

@author: My Nguyen 

Make custom topology for CG simulations

"""
#import mdtraj as md
from simtk import openmm, unit
from simtk.unit import *
from simtk.openmm import app
import numpy as np
import simtk.openmm.app.topology as topology

#========================================
#0) DEFINE ENERGY, LENGTH & MASS SCALES###
#========================================
epsilon  = 1. * kilojoules_per_mole
sigma = 1. * nanometer
mass = 12 * amu
N_av = 6.022140857*10**23 /mole
kb = 1.380649*10**(-23)*joules/kelvin* N_av #joules/kelvin/mol
#====================================
#1) SYSTEM DIMENSIONLESS PARAMETERS###
#====================================
#a)Molecules:
# molecular structures and numbers of copies
mols = [['H']*24 + ['T']*20, ['S']]
mol_names = ["Pol", "Sol"]
nmols = [10, 50000]

charge = 0.0 * elementary_charge
#density in (length)**-3
reduced_density = 6000./(40.**3)

#========================================
#5) Create a system and add particles to it
#========================================
system = openmm.System()
# Particles are added one at a time
# Their indices in the System will correspond with their indices in the Force objects we will add later
for imol, mol in enumerate(mols):
    for i in range(nmols[imol] * len(mol)):
        system.addParticle(mass)
print("Total number of paricles in system: {}".format(system.getNumParticles()))

# Set the periodic box vectors:
number_density = reduced_density / sigma**3
volume = system.getNumParticles() * (number_density ** -1)
box_edge = [volume ** (1. / 3.)] * 3
box_vectors = np.diag([edge/angstrom for edge in box_edge]) * angstroms
system.setDefaultPeriodicBoxVectors(*box_vectors)

#==================
#6) Create topology
#==================
#Topology consists of a set of Chains 
#Each Chain contains a set of Residues, 
#and each Residue contains a set of Atoms.
mols_flat = [item for sublist in mols for item in sublist]
elements = {}
for i,a in enumerate(np.unique(mols_flat)):
    elements.update({a : app.element.Element(200 + i, 'Pol', 'g'+a, mass)})

def makeTop(nmols, mols):
    top = topology.Topology()
    for ispec in range(len(nmols)): #loop over each species
        for imol in range(nmols[ispec]): #loop over each molecule within species
            chain = top.addChain() #create chain for each molecule
            resname = mol_names[ispec]
            res = top.addResidue( resname, chain) # each molecule is one residue
            for atomInd,atomName in enumerate(mols[ispec]):
                el = elements[atomName]
                if atomInd > 0:
                    previousAtom = atom
                atom = top.addAtom( atomName, el, res )
                if atomInd > 0:
                    top.addBond(previousAtom,atom)
    return top
print ("\nCreating topology")
top = makeTop(nmols, mols)

#==============================
#10) Prepare the Simulation ##
#==============================
positions = [box_edge[0]/angstrom,box_edge[1]/angstrom,box_edge[2]/angstrom]*angstrom * np.random.rand(system.getNumParticles(),3)
initialpdb = 'top.pdb'
app.PDBFile.writeModel(top, positions, open(initialpdb,'w'))
