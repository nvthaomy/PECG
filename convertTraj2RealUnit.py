import mdtraj

top = __top__
traj = __traj__
scale = __scale__

t = mdtraj.load(traj,top)
t.unitcell_lengths *= scale
t.xyz *= scale
trajOut = traj.split('.dcd')[0] + 'realUnit.dcd'
t.save(trajOut)
