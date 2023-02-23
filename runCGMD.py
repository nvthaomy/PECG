"""
Run CG MD from Srel forcefield file
"""
from simtk import openmm, unit
from simtk.openmm import app
import simtk.openmm.app.topology as topology
import mdtraj
import time
import numpy as np
import sys, os
sys.path.append('/home/mnguyen/bin/cgfts')
from cgfts.forcefield.forcefield_v2 import ForceField
from collections import OrderedDict
sys.path.append('/home/mnguyen/bin/PEFTS/')
from writePolyFTS import PolyFTS

# force field parameters
ff_file = '../../PDA_ff.dat'
CG_sigma = 0.310741   #nm, set this to 1 if want to use the same length scale as ff file
lb = 0.715231667721/0.310741 * CG_sigma  # To implement in OMM, unit charge is now 0.0848385924569441in CG_sigma unit
kT = 1.
r_kT_kTref = 1.

# bead species
bead_types = ['HOH','Na+', 'Cl-', 'SO4', 'C2', 'C+','COH', 'EO']
charges = [0, 1, -1, -1, 0, 1, 0, 0]

# molecules
chain_names = ['Poly', 'DS', 'C13EO11', 'DS', 'C13EO11', 'Na+', 'Cl-', 'HOH']
nchains = [1, 34, 76, 34, 76, 3536, 3516, 466166]
chains = [['C+']* 48, ['SO4'] +  6*['C2'], ['C2'] * 6 + ['EO'] * 11 + ['COH'], 
['SO4'] +  6*['C2'], ['C2'] * 6 + ['EO'] * 11 + ['COH'], 
['Na+'],['Cl-'],['HOH']] # assume linear structure
#rigid bond
rigid_bonds = {} #{'HOH': [r]}

# box size
Lx, Ly, Lz = [28,28,28] # in CG_sigma unit

# seeds
init_pdb = None
state2load = None #'../sim/PDA0_checkpnt.chk' # .chk or .xml, overwrite init_pdb

# integration parameters
cut_off = 0.63 * 6 # in CG_sigma unit
reduced_dt = 0.01
reduced_temp = 1.
pressure = 8.52 / CG_sigma**3 # kT/nm3
press_ax = None # axis (0,1 or 2) to allow fluctuation in anisotropic barostat

min_steps = 1000
equil_steps = 0 #int(1000/reduced_dt)
steps =  int(50000/reduced_dt)
nframes = 2000
report_freq = int(steps/nframes) 

#Platform settings
platformName = 'CUDA'
platformPrecision = 'mixed' #activated only if device != -1
device = -1 #-1 to let omm automatically choose platform type, device, and precision

#################### End of user inputs ####################
# eliminate bead_types if their occurence is zero:
nbeads = [0] * len(bead_types)
for im,m in enumerate(chains):
    lst = nchains[im] * m
    atoms, counts = np.unique(lst,return_counts=True)
    for ia,a in enumerate(atoms):
        nbeads[bead_types.index(a)] += counts[ia]

im_nonzero = np.where(np.array(nchains) > 0)[0]
if im_nonzero.shape:
    chain_names = [chain_names[im] for im in im_nonzero]
    nchains = [nchains[im] for im in im_nonzero]
    chains = [chains[im] for im in im_nonzero]
ia_nonzero = np.where(np.array(nbeads) > 0) [0]   
if ia_nonzero.shape:
    bead_types = [bead_types[ia] for ia in ia_nonzero]
    charges = [charges[ia] for ia in ia_nonzero]
    nbeads = [nbeads[ia] for ia in ia_nonzero]
nspec = len(bead_types)

# check electroneutrality
qtot = np.sum(np.array(nbeads) * np.array(charges))
if np.abs(qtot) > 1e-3:
    raise Exception('system is not neutral, total charge: {} e'.format(qtot)) 
else:
    print('System is neutral, total charge: {} e'.format(qtot))
  
# import sim ff file
ff = ForceField.from_sim_ff_file(ff_file, kT=kT)
ff.reorder_bead_types(bead_types)
asmear = [bt.smear_length * CG_sigma for bt in ff.bead_types]
print('Bead types from ff file: {}'.format(' '.join([bt.name for bt in ff.bead_types])))
print('smear lengths: {}'.format(asmear))
print('charges: {}'.format(charges))
print('Chain types: {}'.format(chain_names))
print(chains)
print('chain numbers: {}'.format(nchains))
print('total {} atoms'.format(np.sum(nbeads)))
# excluded volume matrix
B = np.zeros((nspec, nspec)) # gaussian prefactor
Kappa = np.zeros((nspec, nspec)) # kappa
# bond
Dist0 = np.zeros((nspec, nspec))
FConst = np.zeros((nspec, nspec))
for i, bt1 in enumerate(bead_types):
    for j in range(i,nspec):
        bt2 = bead_types[j]
        try:
            Dist0[i,j] = Dist0[i,j] = ff.get_pair_potential("Bonded", bt1, bt2).Dist0.value * CG_sigma
            FConst[i,j] = FConst[j,i] = ff.get_pair_potential("Bonded", bt1, bt2).FConst.value * r_kT_kTref / CG_sigma**2
        except:
            pass
        B[i,j] = B[j,i] = ff.get_pair_potential("Gaussian", bt1, bt2).B.value * r_kT_kTref 
        Kappa[i,j] = Kappa[j,i] = ff.get_pair_potential("Gaussian", bt1, bt2).Kappa.value / CG_sigma**2
print('--- Forcefield read from {} ---'.format(ff_file))
s='B:\n{}\n'.format(B)
s+='Kappa:\n{}\n'.format(Kappa)
s+='FConst:\n{}\n'.format(FConst)
s+='Dist0:\n{}\n'.format(Dist0) 
print(s)
with open('ff_check.txt','w') as f:
    f.write(s)

# set up openMM
#========================================
#Some simulation settings
#========================================
constraintTolerance = 1e-6
ewald_tolerance = 1e-5

#========================================
###DEFINE ENERGY, LENGTH & MASS SCALES###
#Try to use OMM's built-in unit definitions to the extent possible#
#========================================
epsilon  = 1.0 * unit.kilojoules_per_mole     #kJ/mole
sigma = 1.0 * unit.nanometer                #nm
tau = 1.0*unit.picoseconds                  #ps
mass = tau**2/sigma**2*epsilon              #dalton = g/mole

N_av = unit.constants.AVOGADRO_CONSTANT_NA #6.02214179e23 / mole
kB = unit.constants.BOLTZMANN_CONSTANT_kB #1.3806504e-23 * joule / kelvin
Tref = epsilon/kB/N_av

#from openmmtools.constants import ONE_4PI_EPS0 
ONE_4PI_EPS0 = 138.935456       #in OMM units, 1/4pi*eps0, [=] kT/mole*nm/e^2
qFactor = ONE_4PI_EPS0**-0.5    #s.t. electrostatic energy (1/4pi eps0)*q^2*qF^2/r = 1kB*Tref*Nav/mole = 1kJ/mole = 1epsilon

print("\n=== OMM Units: ===")
print("mass: {} dalton".format(mass.value_in_unit(unit.dalton)))
print("epsilon: {}".format(epsilon))
print("tau: {}".format(tau))

#========================================
# Create a system and add particles to it
#========================================
print("=== Creating OMM System ===")
system = openmm.System()

# Set the periodic box vectors:
box_edge = [Lx,Ly,Lz]
box_vectors = np.diag(box_edge) * sigma
system.setDefaultPeriodicBoxVectors(*box_vectors)
print('Box vectors:\n{}'.format(system.getDefaultPeriodicBoxVectors()))

#==================
##CREATE TOPOLOGY##
#==================
#Topology consists of a set of Chains 
#Each Chain contains a set of Residues, 
#and each Residue contains a set of Atoms.
#Currently we do 1 residue per chain, and residue is effectively the Molecule class in sim.
print("\n=== Creating Topology ===")
elements = {}
atomTypeIndex = {} # atom name: atom index
atomNameMap = []
for ia,a in enumerate(bead_types):
    newsymbol = 'Z{}'.format(ia)
    if newsymbol not in app.element.Element._elements_by_symbol:
        elements[a]=app.element.Element(200+ia, a, newsymbol, mass)
    else:
        elements[a]=app.element.Element.getBySymbol(newsymbol)
    atomTypeIndex[a] = ia
    atomNameMap.append(a)
print('atom name to atom type index: ',atomTypeIndex)
#print(atomNameMap)

# now create topology
for a in range(np.sum(nbeads)):
    system.addParticle(mass)
print("Total number of particles in system: {}".format(system.getNumParticles()))

top = app.topology.Topology()
top.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
mdtrajtop = app.topology.Topology() #so that later can make molecules whole
# atom and bond list over all molecule copies
atom_list = [] 
bond_list = [] # indices
ia_current = 0
rigid_bond_num = 0
for imt, mol_type in enumerate(chains): #iterate through all chain types 
    mol_name = chain_names[imt]
    mol = chains[imt]
    for im in range(nchains[imt]): #iterate through chain copies
        chain = top.addChain() #Create new chain for each molecule
        res = top.addResidue(mol_name,chain)
        mdt_chain = mdtrajtop.addChain() #Create new chain for each molecule
        mdt_res = mdtrajtop.addResidue(mol_name,mdt_chain)
        
        # add the atoms
        for ia,a in enumerate(mol):
            atom_list.append(a)
            el = elements[a]
            if ia > 0:
                previousAtom = atom
            atom = top.addAtom( a, el, res )
            mdt_atom = mdtrajtop.addAtom( a, mdtraj.element.Element.getByAtomicNumber(atomTypeIndex[a]), mdt_res ) #use a dummy element by matching atomic number == cgAtomTypeIndex
            
            # add the bonds
            if ia > 0: # always bond this atom with the previous
                top.addBond(previousAtom,atom)
                mdtrajtop.addBond(previousAtom,atom) #don't worry about adding constraint to the mdtraj topology

                # add rigid bond, will ignore harmonic bond if detected in forcefield file
                if is_rigid:
                    system.addConstraint( ia_current-1, ia_current, rigid_bonds[mol_name][ia-1])
                    rigid_bond_num += 1
                else: # add harmonic bond later
                    bond_list.append([ia_current-1,ia_current])
            ia_current += 1
print("Total number of constraints: {}".format(system.getNumConstraints()))   
print('{} rigid bonds'.format(rigid_bond_num)) 
#convert atom_list to dictionary of atom indices for each bead species
atom_dict = {aname: np.squeeze(np.where(np.array(atom_list) == aname)) for aname in bead_types}

#====================
##CREATE FORCEFIELD##
#====================
# currently allow: Gaussian interaction, harmonic bonds, smear Coulomb
print("\n=== Create Forcefield ===")
print("---> Harmonic bond")
bondedForce = openmm.HarmonicBondForce()
for ib,(i,j) in enumerate(bond_list): 
    ai = atomTypeIndex[atom_list[i]] # index in bead_types
    aj = atomTypeIndex[atom_list[j]]
    dist0 = Dist0[ai,aj]
    fconst = FConst[ai,aj]
    bondedForce.addBond(i, j, dist0*sigma, 2.0*fconst*epsilon/sigma/sigma)
print('{} bonds'.format(len(bond_list)))
system.addForce(bondedForce)

print("---> Nonbonded Force") #Gaussian + ShortRange (Including ewald smearing correction)+ ExtField#
nonbondedMethod = 2
print("    -> Gaussian")
B *= epsilon.value_in_unit(unit.kilojoule/unit.mole)
Kappa = Kappa # already in nm^-2
energy_function = 'B(type1,type2)*exp(-Kappa(type1,type2)*r^2);'
fcnb = openmm.CustomNonbondedForce(energy_function)
fcnb.addPerParticleParameter('type')
fcnb.setCutoffDistance(cut_off)
fcnb.setNonbondedMethod( nonbondedMethod ) #2 is cutoff non periodic
fcnb.addTabulatedFunction('B', openmm.Discrete2DFunction(nspec,nspec,B.ravel(order='F')) )
fcnb.addTabulatedFunction('Kappa', openmm.Discrete2DFunction(nspec,nspec,Kappa.ravel(order='F')) )
for atom in top.atoms():
    fcnb.addParticle( [atomTypeIndex[atom.name]] )
system.addForce(fcnb)

if  any(charges):
    print('    -> Ewald')
    nbfmethod = openmm.NonbondedForce.PME 
    print("    To implement in OMM, unit charge is now {}".format(qFactor))
    chargeScale = qFactor * lb **0.5

    nbf = openmm.NonbondedForce()
    nbf.setCutoffDistance(cut_off)
    nbf.setEwaldErrorTolerance( ewald_tolerance )
    nbf.setNonbondedMethod( nbfmethod )
    nbf.setUseDispersionCorrection(False)
    nbf.setUseSwitchingFunction(False)

    for (i,atom) in enumerate(atom_list):
        charge = charges[atomTypeIndex[atom]] * chargeScale #In dimensionless, EwaldPotential.Coef is typically 1, and usually change relative strength via temperature. But can also scale the coef, which then acts as lB in the unit length
        LJsigma = 1.0
        LJepsilon = 0.0
        nbf.addParticle(charge, LJsigma, LJepsilon)
    system.addForce(nbf)

    print('    -> Ewald smear correction')
    aborn_matrix = np.zeros((nspec,nspec))  # = asmear * sqrt(pi) = sqrt((a1**2 + a2**2)/2) * sqrt(pi)
    for i,sp1 in enumerate(bead_types):
        for j,sp2 in enumerate(bead_types):
            aborn_matrix[i,j] = np.sqrt((asmear[i]**2. + asmear[j]**2.)/2.) * np.sqrt(np.pi)
    energy_function = 'coef*q1*q2 * ( (erf(factor*r) - 1)/r - shift );'
    energy_function += 'shift = (erf(factor*rcut) -1)/rcut;'
    energy_function += 'factor = sqrt({:f})/2/aborn(type1,type2);'.format(np.pi)
    energy_function += 'coef = {:f};'.format(lb)
    energy_function += 'rcut = {:f};'.format(cut_off)
    fcnb = openmm.CustomNonbondedForce(energy_function)
    fcnb.addPerParticleParameter('type')
    fcnb.addPerParticleParameter('q')
    fcnb.setCutoffDistance( cut_off )
    fcnb.setNonbondedMethod( openmm.NonbondedForce.CutoffPeriodic )
    fcnb.addTabulatedFunction('aborn',openmm.Discrete2DFunction(nspec, nspec, aborn_matrix.ravel(order='F')))
    for (i,atom) in enumerate(atom_list):
        q = charges[atomTypeIndex[atom]]
        fcnb.addParticle( [atomTypeIndex[atom], q] )
    system.addForce(fcnb)

print("--- All forces added: ---")
forces = [f for f in system.getForces()]
for f in forces:
    print(f)

#add box vector in top to make sure box info is written
top.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
print('Uses periodic: {}'.format( system.usesPeriodicBoundaryConditions() ))

#=======================
## PREPARE SIMULATION ##
#=======================
dt = reduced_dt * tau
temperature = reduced_temp * epsilon/kB/N_av
friction = 1/(100 * dt)

if pressure:
    useNPT = True
    useNVT = False
else:
    useNPT = False
    useNVT = True
"""
reduced_pressure = 1
reduced_Pdamp = 0.1 #time units
pressure = reduced_pressure * epsilon/(sigma**3) / N_av
barostatInterval = int(reduced_Pdamp/reduced_dt)
"""
print("\n=== Preparing Simulation ===")
if useNPT:
    pressure = pressure * epsilon/N_av/sigma/sigma/sigma #convert from unitless to OMM units
    barostatInterval = 25 #in units of time steps. 25 is OpenMM default

    if press_ax is None:
        system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
        print("Isotropic Barostat")
    elif press_ax in [0,1,2]:
        Pscale = [False,False,False]
        Pscale[ press_ax] = True
        system.addForce(openmm.MonteCarloAnisotropicBarostat(openmm.vec3.Vec3(pressure,pressure,pressure), temperature, Pscale[0], Pscale[1], Pscale[2], barostatInterval))   
        print("Anisotropic Barostat, axes fluctuations: {}".format(Pscale))
    print("Added MC Barostat with P {} (eps/sig^3), T {}, freq {}".format(
        6.02214179e23 * pressure.value_in_unit(unit.kilojoules/unit.nanometer**3),temperature,barostatInterval))
    print("In OMM units, is P={}'bar'".format(pressure.value_in_unit(unit.bar)))

integrator = openmm.LangevinIntegrator(temperature, friction, dt)
if system.getNumConstraints() > 0:
    print("Applying bond constraints before starting")
    integrator.setConstraintTolerance(constraintTolerance)
if device == -1:
    print("Automatically choosing platform")
    simulation = app.Simulation(top,system,integrator)
else:
    print("Manually setting platform: {} and device: {}".format(platformName,device))
    platform = openmm.Platform.getPlatformByName(platformName)
    if platformName == 'CUDA':
        platform.setPropertyDefaultValue("CudaDeviceIndex", str(device))
        platform.setPropertyDefaultValue("Precision",str(platformPrecision))
    elif platformName == 'OpenCL':
        platform.setPropertyDefaultValue("OpenCLDeviceIndex", str(device))
        platform.setPropertyDefaultValue("Precision",str(platformPrecision))
    else:
        platform = openmm.Platform.getPlatformByName('CPU')
    simulation = app.Simulation(top,system, integrator, platform)
    print('Default precision: {}'.format(platform.getPropertyValue(simulation.context,'Precision')))

simOptions = {'dt':dt, 'temp':temperature, 'fric':friction}
print(simOptions)

# prepare log files
chkfile = 'checkpnt.chk'
simulation.reporters.append(app.checkpointreporter.CheckpointReporter(chkfile,10000))
trajfile = 'traj.dcd'
eq_logfile = "eqlog.txt"
logfile = "prodlog.txt"

print("\n=== Initializing Simulation ===")
if init_pdb: 
    pos = app.PDBFile(init_pdb).getPositions(asNumpy = True)
else:
    pos = [Lx,Ly,Lz] * np.random.rand(np.sum(nbeads),3) # same unit as Lx, Ly, Lz
pos = pos[np.newaxis,:]
print('positions shape ',np.array(pos).shape)
#need to make molecules whole
unitcell_lengths = np.array([[Lx,Ly,Lz]])
unitcell_angles = np.array([[90., 90., 90.]])

mdt_traj = mdtraj.Trajectory(pos, topology = mdtrajtop, unitcell_lengths = unitcell_lengths, unitcell_angles = unitcell_angles)
bondlist = np.array([ [b[0].index,b[1].index] for b in top.bonds() ], dtype=np.int32) #should be sorted, as long as the bonds come out sorted from the sim Sys object
if len(bondlist) > 0:
    print(">0 bonds, making trajectory whole to be safe")
    wholetraj = mdt_traj.make_molecules_whole(inplace=False,sorted_bonds=bondlist)
    simulation.context.setPositions(wholetraj.xyz[0])#Pos)
else:
    simulation.context.setPositions(pos[0])
initstatefile =  'output.xml'
simulation.saveState(initstatefile)
initialpdb = "initial.pdb"
initial_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
app.PDBFile.writeModel(simulation.topology, initial_positions, open(initialpdb,'w'))

if system.getNumConstraints() > 0:
    simulation.context.applyConstraints(constraintTolerance)
    constrainedpdb = "constrained.pdb"
    constrained_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    app.PDBFile.writeModel(simulation.topology, constrained_positions, open(constrainedpdb,'w'))

#=====================
## RUN SIMULATION ##
#=====================
print('\n=== Run simulations ===')
if state2load:
    try:
        simulation.loadState(state2load)
        print('Initiate simulation from {}'.format(state2load))
    except:
        simulation.loadCheckpoint(state2load)
        print('Initiate simulation from {}'.format(state2load))
    positions = simulation.context.getState(getPositions=True).getPositions()
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    box = simulation.context.getState().getPeriodicBoxVectors()
    print('Current box {}'.format(box))
else:
    if min_steps > 0:
        print(" ... Running energy minimization ...")
        tmp_pos = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
        app.PDBFile.writeModel(simulation.topology, tmp_pos, open("pre-minimized.pdb",'w'))

        minimizefile = "minimized.pdb"
        simulation.minimizeEnergy(maxIterations=min_steps)
        simulation.context.setVelocitiesToTemperature(simOptions["temp"]*3)
        minimized_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
        app.PDBFile.writeModel(simulation.topology, minimized_positions, open(minimizefile,'w'))

app.PDBFile.writeModel(simulation.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),open('top.pdb','w'))
app.PDBFile.writeFooter(simulation.topology,open('top.pdb','a'))

print(' ... Equilibrating ...')
simulation.reporters.append(app.StateDataReporter(sys.stdout, report_freq*100, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
simulation.reporters.append(app.StateDataReporter(eq_logfile, report_freq, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
print("Progress will be reported every {} steps".format(report_freq))
simulation.step(equil_steps)
equilibratefile = "equilibrated.pdb"
equilibrated_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
app.PDBFile.writeModel(simulation.topology, equilibrated_positions, open(equilibratefile,'w'))

print(' ... Production ...')
print('{} steps total, {} tau total ({})'.format(steps, steps * reduced_dt, steps * reduced_dt *tau ))
simulation.reporters.pop()
simulation.reporters.append(app.StateDataReporter(logfile, report_freq, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
simulation.reporters.append(mdtraj.reporters.DCDReporter(trajfile,report_freq))

blocksteps = min(steps,10000)
nblock = int(steps/blocksteps)
start = time.time()
for i in range(nblock):
    simulation.step(blocksteps)
    simulation.saveState('output.xml')
end = time.time()
equilibratefile = "final.pdb"
equilibrated_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
app.PDBFile.writeModel(simulation.topology, equilibrated_positions, open(equilibratefile,'w'))
print("=== Finished production in {} seconds".format(end-start))
