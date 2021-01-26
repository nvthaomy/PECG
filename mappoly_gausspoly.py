#!/usr/bin/env python
import sys, pickleTraj,time, os
import numpy as np
import mdtraj as md
def getMap(traj, top, nameMap, mapratio=1):
    """map AA traj to CG traj
    assuming map one AA residues to one CG bead
    input: nameMap = {AAres:CGatomType}
    outputs: [CGatomTypes], [[AtomIds in bead 1],[AtomIds in bead 2], ...]"""    
    Mass1List = []
    AAatomId = []
    AAres = []
    CGatomTypes = []
    CGchainId = [] #chain id for each corresponding CGatomTypes
    for res in top.residues:
        AAres.append(res.name)
        if res.name == 'Pol':
            if float(res.index) % float(mapratio) == 0:
                CGatomTypes.append(nameMap[res.name])
                CGchainId.append(top.residue(res.index).chain.index)
                #get atom indices of all atoms in this residue
                Atoms1 = []
                Mass1 = []
            for atom in res.atoms:
                Atoms1.append(atom.index)
                Mass1.append(1) #atom.element.mass)
            if float(res.index) % float(mapratio) == 0:
                AAatomId.append(Atoms1)
                Mass1List.append(Mass1)

        else: # always 1:1 mapping for solvent
            CGatomTypes.append(nameMap[res.name])
            CGchainId.append(top.residue(res.index).chain.index)
            Atoms1 = []
            Mass1 = []
            for atom in res.atoms:
                Atoms1.append(atom.index)
                Mass1.append(1) #atom.element.mass)
            AAatomId.append(Atoms1)
            Mass1List.append(Mass1)
    Mol = []
    for bond in top.bonds:
        atom1,atom2 = bond[0],bond[1]
        res1,res2 = atom1.residue, atom2.residue
        if not res1 in Mol:
            Mol.append(res1)
        if res1 != res2 and not res2 in Mol:
            Mol.append(res2)
    return AAatomId, CGatomTypes, AAres, Mass1List,CGchainId 

def convertTraj(traj, top, lengthScale = 1., stride = 1, outTrajExt = '.lammpstrj'):
    """convert trajectory to a specified format and scale box with lengthScale"""
    import mdtraj as md
    import os,shutil
    #to make traj in the current dir
    cwd = os.getcwd()
    outTraj = traj.split('/')[-1]
    outTraj = '.'.join(outTraj.split('.')[:-1]) + outTrajExt
    traj = md.load(traj, top = top, stride = stride)
    if lengthScale != 1.:
        print('Scaling positions and box size by 1/{}'.format(lengthScale))
    traj.xyz /= lengthScale
    traj.unitcell_lengths /= lengthScale
    traj.save(outTraj)
    shutil.move(outTraj,os.path.join(cwd,outTraj))
    print('moving traj to {}'.format(os.path.join(cwd,outTraj))) 
    outTraj = os.path.join(cwd,outTraj)
    return outTraj

def mapTraj(traj, top, nameMap, lengthScale, stride=1, mapratio=1, wrapTraj=False):
    import sim
    """ Assume polymers are listed first
        CGatomTypes: 1 by n list of CG atom types
        AAatomID: n by x list of indices of AA atoms in CG beads
        lengthScale (nanometer): divide positions and box dimensions by this scale, for conversion between real and dimensionless units"""

    traj1 = md.load(traj, top = top, stride=stride)
    _top1 = traj1.topology

    AAatomId, CGatomTypes, AAres, Mass1List, CGchainId= getMap(traj1, _top1, nameMap, mapratio=mapratio )
    AtomTypes, counts = np.unique(CGatomTypes, return_counts = True)

    print("\n ===== CG Atom count ===== :")
    for i, atom in enumerate(AtomTypes):
        print('%s\t%i'%(atom,counts[i]))
    
    if lengthScale != 1.:
        print('Scaling positions and box size by 1/{}'.format(lengthScale))
    traj1.xyz /= lengthScale
    traj1.unitcell_lengths /= lengthScale
    
    outTraj = traj.split('/')[-1]
    outTraj = '.'.join(outTraj.split('.')[:-1])
    
    # setup topology file for COM trajectory
    element = md.element.Element(200,"gaussPolymer","gP", 1., 1.)
    _top2 = md.Topology()
    
    for cindx in range(_top1.n_chains):
         _top2.add_chain()
         rindx_in_chain = np.where(np.abs(np.array(CGchainId,dtype=float)-float(cindx))<1e-4)[0]
         for ii,resname in enumerate(np.array(CGatomTypes)[rindx_in_chain]):
             _top2.add_residue(resname,_top2.chain(-1))
             _top2.add_atom(resname,element,_top2.residue(-1))
             if ii > 0:
                 _top2.add_bond(_top2.atom(-2),_top2.atom(-1))

    print("\n ===== Mapping and Writing Trajectory =====")
    if wrapTraj:
        print('Wrap trajectory')
    start = time.time()
    # setup np matrix to hold the mapped coordinates
    COM_xyz = np.zeros((int(traj1.n_frames),int(len(CGatomTypes)),3))
    # loop through each residue and do mapping for ALL frames in trajectory
    CGrindx = 0
    for rindx, residue in enumerate(_top1.residues):
        if residue.name == 'Pol':
            if float(rindx) % float(mapratio) == 0:
                resTraj = traj1.atom_slice(_top1.select('resid {} to {}'.format(rindx,rindx+mapratio-1)))
                tempmass = np.column_stack((1.,1.,1.))
                rsq = np.multiply(resTraj.xyz,tempmass)
                rsq = np.sum(rsq,axis=1) / float(resTraj.n_atoms)
                if wrapTraj: # WrapTrajCoordinates
                    rsq = np.mod(rsq,traj1.unitcell_lengths)
                COM_xyz[:,CGrindx,:] = rsq        
                CGrindx+=1
        else:
            resTraj = traj1.atom_slice(_top1.select('resid {}'.format(rindx)))
            tempmass = np.column_stack((1.,1.,1.))
            rsq = np.multiply(resTraj.xyz,tempmass)
            rsq = np.sum(rsq,axis=1) / float(resTraj.n_atoms)
            if wrapTraj: # WrapTrajCoordinates
                rsq = np.mod(rsq,traj1.unitcell_lengths)
            COM_xyz[:,CGrindx,:] = rsq
            CGrindx+=1
    final = time.time()
    totaltime = final - start
    print('RunTime:            {}\n'.format(totaltime))
    
    MappedTrj = md.Trajectory(COM_xyz,_top2,unitcell_lengths = traj1.unitcell_lengths, unitcell_angles = traj1.unitcell_angles)
    MappedTrj.save_lammpstrj(os.path.join(os.getcwd(),outTraj+"_mapped.lammpstrj"))
    MappedTrj.save_dcd(os.path.join(os.getcwd(),outTraj+"_mapped.dcd"))
    _PDB = md.formats.PDBTrajectoryFile(os.path.join(os.getcwd(),outTraj+"_mapped.pdb"), mode='w', force_overwrite=True, standard_names=True)
    _PDB.write(MappedTrj.xyz[0],MappedTrj.topology)    
 
    # Convert to sim format
    MappedTrj = pickleTraj(outTraj+"_mapped.lammpstrj")
    BoxL = MappedTrj.FrameData['BoxL']
    print('Box L {}'.format(BoxL))
    return  CGatomTypes, MappedTrj, BoxL        
    
if __name__ == "__main__":
    import argparse as ap
    parser = ap.ArgumentParser(description="mapping AA trajectory")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    parser.add_argument('-ls', type = float, help = "scale traj by 1/ls", default =1)
    parser.add_argument('-m', default=1,type=int, help = 'AA bead per CG bead, mapping ratio')
    parser.add_argument('-w', action='store_true', help = 'Wrap trajectory')
    args = parser.parse_args()
    
    traj = sys.argv[1]
    top = sys.argv[2]
    stride = args.stride
    lengthScale = args.ls
    mapratio = args.m
    wrapTraj = args.w

    #dictionary defines: {AAres:[CGatomTypes]}
    nameMap = {'Pol':'A', 'Sol':'B'}
    mappedTraj = mapTraj(traj, top, nameMap, lengthScale, stride = stride, mapratio = mapratio, wrapTraj=wrapTraj)

