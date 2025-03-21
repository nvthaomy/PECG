#!/usr/bin/env python
import sys, pickleTraj
import numpy as np
def getMap(traj, top, nameMap):
    """map AA traj to CG traj
    assuming map one AA residues to one CG bead
    input: nameMap = {AAres:CGatomType}
    outputs: [CGatomTypes], [[AtomIds in bead 1],[AtomIds in bead 2], ...]"""    
    import mdtraj as md
    traj = md.load(traj, top = top)
    top = traj.topology
    Mass1List = []
    AAatomId = []
    AAres = []
    CGatomTypes = []
    for res in top.residues:
        if nameMap[res.name] == 'skip': 
            print('Skip mapping residue {}'.format(res.name))
            continue
        AAres.append(res.name)
        CGatomTypes.append(nameMap[res.name])
        #get atom indices of all atoms in this residue
        Atoms1 = []
        Mass1 = []
        for atom in res.atoms:
            Atoms1.append(atom.index)
            Mass1.append(atom.element.mass)
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
            
    return AAatomId, CGatomTypes, AAres, Mass1List 

def convertTraj(traj, top, lengthScale = 1., stride = 1, start_frame=0, outTrajExt = '.lammpstrj'):
    """convert trajectory to a specified format and scale box with lengthScale"""
    import mdtraj as md
    import os,shutil
    #to make traj in the current dir
    cwd = os.getcwd()
    outTraj = traj.split('/')[-1]
    outTraj = '.'.join(outTraj.split('.')[:-1]) + outTrajExt
    traj0 = md.load(traj, top = top)
    traj = traj0[start_frame:]
    traj = traj[::stride]
    if lengthScale != 1.:
        print('Scaling positions and box size by 1/{}'.format(lengthScale))
    traj.xyz /= lengthScale
    traj.unitcell_lengths /= lengthScale
    traj.save(outTraj)
    shutil.move(outTraj,os.path.join(cwd,outTraj))
    print('moving traj to {}'.format(os.path.join(cwd,outTraj))) 
    outTraj = os.path.join(cwd,outTraj)
    return outTraj

def mapTraj(traj, top, nameMap, lengthScale, stride=1, start_frame=0, outTrajExt='.lammpstrj'):
    import sim
    """ CGatomTypes: 1 by n list of CG atom types
        AAatomID: n by x list of indices of AA atoms in CG beads
        lengthScale (nanometer): divide positions and box dimensions by this scale, for conversion between real and dimensionless units"""

    assert outTrajExt in ['.lammpstrj','.pdb'] # only allow these extensions

    AAatomId, CGatomTypes, AAres, Mass1List = getMap(traj, top, nameMap)
    AtomTypes, counts = np.unique(CGatomTypes, return_counts = True)
    # ===== create mapped object =====
    print("\n ===== Creating Index Mapper =====")
    Map = sim.atommap.PosMap()
    for i, CGatomType in enumerate(CGatomTypes):
        Atoms1 = AAatomId[i]
        Atom2 = i
        Mass1 = Mass1List[i]
        this_Map = sim.atommap.AtomMap(Atoms1 = Atoms1, Atom2 = Atom2, Mass1=Mass1)
        Map += [this_Map]
    
    print("\n ===== Converting AA traj to lammpstrj format  =====")
    traj = convertTraj(traj, top, lengthScale = lengthScale, stride = stride, start_frame=start_frame, outTrajExt='.lammpstrj')
    if outTrajExt=='.lammpstrj':
        outTraj = traj.split('.lammpstrj')[0] + '_mapped.lammpstrj.gz'    
    elif outTrajExt=='.pdb':
        outTraj = traj.split('.lammpstrj')[0] + '_mapped.pdb'
    # ===== read AA traj =====
    print("\n ===== Reading scaled AA Traj =====")
    Trj = pickleTraj(traj)
    BoxL = Trj.FrameData['BoxL']
    
    # ===== write out new mapped traj =====
    print("\n ===== Mapping and Writing Trajectory =====")
    MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = CGatomTypes, BoxL = BoxL)
    #MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes, BoxL = BoxL)
    print('mapped traj name {}'.format(outTraj)) 
    # ===== convert to lammps =====
    print("\n ===== Converting to LAMMPS =====")
    if outTrajExt=='.lammpstrj':
        sim.traj.base.Convert(MappedTrj, sim.traj.LammpsWrite, FileName = outTraj, Verbose = True)    
    elif outTrajExt=='.pdb':
        sim.traj.base.Convert(MappedTrj, sim.traj.PdbWrite, FileName = outTraj, Verbose = True)
    print("\nCG Atom\tcount:")
    for i, atom in enumerate(AtomTypes):
        print('%s\t%i'%(atom,counts[i]))
    return  CGatomTypes, AAatomId, MappedTrj, BoxL 
    
if __name__ == "__main__":
    import argparse as ap
    parser = ap.ArgumentParser(description="mapping AA trajectory")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    parser.add_argument('-ls', type = float, help = "scale traj by 1/ls", default =1)
    args = parser.parse_args()
    
    traj = sys.argv[1]
    top = sys.argv[2]
    stride = args.stride
    lengthScale = args.ls

    #dictionary defines: {AAres:[CGatomTypes]}
    nameMap = {'Na+':'Na+', 'Cl-':'Cl-', 'HOH': 'HOH', 'WAT': 'HOH', 
               'ATP':'A', 'AHP':'A', 'AP': 'A', 'ATD': 'A-', 'AHD': 'A-', 'AD': 'A-',
               'NTP':'B+', 'NHP':'B+', 'NP': 'B+', 'NTD': 'B', 'NHD': 'B', 'ND': 'B'}
    mappedTraj = mapTraj(traj, top, nameMap, lengthScale, stride = stride)
