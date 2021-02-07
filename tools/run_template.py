# Script to read in single-molecule structures, prepare pdb box, and relevant topology output
# Two operation modes:
#   1) write out a full pdb with connectivity specific to this system
#   2) write out an easily-modifiable .txt file that this function can then use to generate an openmm topology object
#
# example usage:
#   python run.py -L 10.0 --fflist ./TrappeUA_Styrene_Gromos.xml
#   python run.py -prefix run1 -sysxml run0_system.xml -chkstate run0_checkpoint.xml 
#   python run.py -prefix run2 -sysxml run0_system.xml -chkstate run0_checkpoint.chk
#   python run.py -prefix run2 -sysxml run0_system.xml -chkstate run0_checkpoint.chk -PME -tail -top packing_topology.p
#   python run.py -prefix run2 -sysxml run0_system.xml -chkstate run0_checkpoint.chk -PME -tail -top packing_topology.csv packing_topology_bond.dat
#   python run.py -prefix run2 -sysxml run0_system.xml -chkstate run0_checkpoint.chk -PME -tail -top topology.txt
#   python run.py -initpdb packing_post_minimization.pdb -boxfile packing_box.txt -ff ~/mylib/SCOUTff/TrappeUA/XML/TrappeUA_Master.xml -PME -tail -top topology.txt
#
# key variables:
#   sys_descrip: [ ('mol1filename',num_mol1), ('mol2filename',num_mol2), ... ]
#   sys: [ (mol1_ommtopology,num_mol1), (mol2_ommtopology,num_mol2), ... ]
#   top: the openmm topology with desired # molecules
#
# key outputs:
#   topology.txt
#   proposed_initial_packing.pdb
#

from simtk import unit
from simtk.openmm import app
import simtk.openmm as mm
import numpy as np
import argparse as ap
import time
import os
from subprocess import call
import pickle, pandas, mdtraj

def write_full_pdb(filename,topology,positions):
    app.pdbfile.PDBFile.writeHeader(topology,open(filename,'w'))
    app.pdbfile.PDBFile.writeModel(topology,positions,open(filename,'a'))
    app.pdbfile.PDBFile.writeFooter(topology,open(filename,'a'))

# === Simulation Options ===
pressure        = 1.0  # bar
temperature     = __TEMP__ # kelvin
temperature_anneal = __TEMPANNEAL__ #kelvin
barostatfreq    = 25
use_gpu         = True

#nonbonded_method = app.LJPME #set below
nonbonded_cutoff = 1.*unit.nanometer
ewald_tol = 0.0001
friction = 1./unit.picosecond
dt = 0.002 #unit.picosecond

dcd_report_freq = int(__DCDFREQ__./dt)
thermo_report_freq = int(__THERMOFREQ__/dt)

#Annealing
annealing_steps = int(__ANNEALTIME__/dt)
equilibration_steps = int(__EQTIME__/dt)
production_steps = int(__PRODUCTIONTIME__/dt)
_C1 = app.element.Element(201, "C1", "CC", 13.018*unit.dalton)
_C2 = app.element.Element(202, "C2", "CT", 13.018*unit.dalton)


# === Parse ===
parser = ap.ArgumentParser(description="Generate a simulation's topology file from single-molecule pdbs")
parser.add_argument('-prefix', type=str, default='run0', help = 'prefix to collect files')
parser.add_argument('-initpdb', type=str, default='packing_proposal.pdb', help = "initial config. MUST have connectivity.")
parser.add_argument('-L', type=float, default=None, help = '(cubic) box L in nm')
parser.add_argument('-boxfile', type=str, default='packing_box.txt')
parser.add_argument('-ff','--fflist', type=str, default=None, nargs='+', help = "string of path to ff files")
parser.add_argument('-sysxml', type=str, default=None, help = 'system xml')
parser.add_argument('-chkstate', type=str, default=None, help = 'checkpoint or state to load')
parser.add_argument('-PME', action='store_true', help = 'LJPME is default, toggle to enable PME')
parser.add_argument('-tail', action='store_true', help = 'tail correction is false by default, toggle to turn on')
parser.add_argument('-NVT', action='store_true', help = 'NVT, turn off barostat')
parser.add_argument('-top', default=None, type=str, nargs='+', help = 'topology files')
args = parser.parse_args()

run_npt         = not args.NVT
prefix = args.prefix
if args.L is not None:
    print('CAUTION: (cubic) box L {} entered, overriding boxfile {}'.format(args.L, args.boxfile))
    box_L = args.L
    periodic_box_vectors = [[box_L,0.,0.],[0.,box_L,0.],[0.,0.,box_L]]
else:
    periodic_box_vectors = np.loadtxt( args.boxfile ) 

ff_list = args.fflist
pdb = app.PDBFile 

print('args.top: {}'.format(args.top))
if args.top is None:
    print('reading topology from args.initpdb {}'.format(args.top))
    sys_pdb = app.PDBFile(args.initpdb)
    top = sys_pdb.topology
else: #is a list
    print('reading topology from args.top {}'.format(args.top))
    suffixes = [ name.split('.')[-1] for name in args.top ]
    if suffixes[0] == 'p':
        print('... reading a pickle file (should be openmm topology object)')
        top = pickle.load(open(args.top[0],"rb"))
    elif suffixes[0] == 'txt':
        print('... reading our custom topology definition. needs a topology library, construct manually')

        structlib='/home/dans/SCOUTff/tools/../structlib/Kraton/'
        def pdb_to_xyz(filename):
            _t = mdtraj.load(filename)
            prefix,suffix = os.path.splitext(filename)
            newfilename = prefix+'.xyz'
            _t.save(newfilename)
            return newfilename


        with open('topology.txt','r') as f:
            sys_descrip = []
            for line in f:
                tmp = line.split()
                sys_descrip.append( [tmp[0], int(tmp[1])] )
        print(sys_descrip)

        sys = []
        nMols = []
        chainPDBs = []
        for entry in sys_descrip:
            print(entry)
            filetype = entry[0].split('.')[-1]
            pdb = app.PDBFile('{}{}'.format(structlib,entry[0]))
            sys.append( (pdb.topology, entry[1]) )
            if filetype == 'xyz':
                tmp_pdbfile = os.path.join(structlib, entry[0])
                chainPDBs.append( pdb_to_xyz(tmp_pdbfile) )
            else:
                chainPDBs.append(os.path.join(structlib,entry[0]))
            nMols.append(entry[1])

        top = app.topology.Topology()
        for ii,mol in enumerate(sys):
            moltop, nummol = mol
            for jj in range(nummol):
                atoms_in_top = []
                for c in moltop.chains():
                    chain = top.addChain()
                    for r in c.residues():
                        residue = top.addResidue(r.name,chain)
                        for a in enumerate(r.atoms()):
                            atom = top.addAtom(a[1].name, a[1].element, residue)
                            atoms_in_top.append(atom)
                            #print(atom)
                for bond in moltop.bonds():
                    bi1,bi2 = bond[0].index, bond[1].index
                    top.addBond( atoms_in_top[bi1], atoms_in_top[bi2] ) 
    elif len(suffixes) == 2:
        print('... first file is data frame, second file can be read as nparray of bonds')
        df = pandas.read_csv(args.top[0])
        bonds = np.loadtxt(args.top[1])
        top_mdtraj = mdtraj.Topology.from_dataframe(df,bonds)
        top = top_mdtraj.to_openmm()
    else:
        raise ValueError('Unrecognized topology input file formats')

if ff_list is None and args.sysxml is None:
    raise ValueError('must either provide forcefield or system xml')

if args.PME:
    print('using PME with tail correction {}'.format(args.tail))
    nonbonded_method = app.PME
    use_tail = args.tail
    print('Actually, for noncharged PME, can just use cutoff periodic')
    nonbonded_method = app.CutoffPeriodic
else: 
    print('using LJPME, tail correction automatically off')
    nonbonded_method = app.LJPME
    use_tail = False


# === Set up simulation ===
if ff_list is not None:
    forcefield = app.ForceField(*ff_list)
    unmatched_residues = forcefield.getUnmatchedResidues(top)
    print('\n=== Unmatched residues ===\n')
    print(unmatched_residues)


print('\n=== Periodic Box ===')
print(periodic_box_vectors)
top.setPeriodicBoxVectors(periodic_box_vectors*unit.nanometer)


print('\n=== Making System ===')
if args.sysxml is not None:
    print('--- Loading {} ---'.format(args.sysxml))
    system = args.sysxml
else: #ff_list should be defined
    print('--- Making new system with given ff and settings ---')
    system = forcefield.createSystem(top, 
                                    nonbondedMethod = nonbonded_method,
                                    nonbondedCutoff = nonbonded_cutoff, 
                                    ewaldErrorTolerance=ewald_tol, 
                                    rigidWater=True, 
                                    constraints=None)
                                    #constraints=app.AllBonds)

    barostat = mm.MonteCarloBarostat( pressure*unit.bar, temperature_anneal*unit.kelvin, barostatfreq )
    #barostat = mm.MonteCarloBarostat( pressure, temperature, barostatfreq )
    if run_npt:
        print('NPT run')
        system.addForce(barostat)
    else:
        print('NVT run')


integrator = mm.LangevinIntegrator(temperature_anneal*unit.kelvin, friction, dt*unit.picosecond)
if use_gpu:
    platform = mm.Platform.getPlatformByName('OpenCL')
    properties = {'DeviceIndex':'0', 'Precision':'mixed'}
else:
    platform = mm.Platform.getPlatformByName('CPU')
    properties = {'Threads': '1'}

simulation = app.Simulation( top, system, integrator, platform, properties )

# === Get initial configuration === 
if args.chkstate is not None:
    print("setting system state to {}, but still continuing with prescribed minimization, eq, production protocol".format(args.chkstate))
    suffix = args.chkstate.split('.')[-1]
    if suffix == 'chk':
        simulation.loadCheckpoint(args.chkstate)
    elif suffix == 'xml':
        simulation.loadState(args.chkstate)
    else:
        raise ValueError('unsupported file format .{}'.format(suffix))
else:
    print("initializing system to {}".format(args.initpdb))
    suffix = args.initpdb.split('.')[-1] 
    if suffix == 'pdb':
        sys_pdb = app.PDBFile(args.initpdb)
        positions = sys_pdb.positions
    elif suffix == 'xyz':
        data = np.loadtxt(args.initpdb, skiprows=2, usecols=(1,2,3))
        positions=data
    else:
        raise ValueError('unsupported coordinate file format {}'.format(args.initpdb))
    
    simulation.context.setPositions(positions)
    simulation.context.setPeriodicBoxVectors(periodic_box_vectors[0],periodic_box_vectors[1],periodic_box_vectors[2]) # Set the periodic box vectors
    #app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_initial.pdb'.format(prefix),'w'))
    simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
    write_full_pdb('{}_initial.pdb'.format(prefix),simulation.topology,positions)

    #simulation.context.applyConstraints(1e-8)


# By default PME turns on tail correction. Manually turn off if requested
print('system may deviate from loaded system, because setting tail correction to {}'.format(use_tail))
system = simulation.context.getSystem() 
ftmp = [f for ii, f in enumerate(system.getForces()) if isinstance(f,mm.NonbondedForce)]
fnb = ftmp[0]
fnb.setUseDispersionCorrection(use_tail)

# write out xml, ONLY after making sure no more changes are made to the system
from simtk.openmm import XmlSerializer
serialized_system_gromacs = XmlSerializer.serialize(system)
outfile = open('{}_system.xml'.format(prefix),'w')
outfile.write(serialized_system_gromacs)
outfile.close()


# === Minimize Energy ===
print('\n=== Minimizing ===')
time_start = time.time()
print('pre minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy(tolerance=0.01*unit.kilojoules_per_mole,maxIterations=10000)
print('post minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
time_end = time.time()
print("done with minimization in {} minutes\n".format((time_end-time_start)/60.))
positions = simulation.context.getState(getPositions=True).getPositions()       
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_minimization.pdb'.format(prefix),'w'))
simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
write_full_pdb('{}_post_minimization.pdb'.format(prefix),simulation.topology,positions)


# === Setup Runs ===

# === Annealing ===
print('\n=== Annealing ===')
simulation.reporters.append(app.statedatareporter.StateDataReporter('{}_thermo_annealing.out'.format(prefix), thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
print('\nannealing at temperature to {} '.format(temperature_anneal))
time_start = time.time()
simulation.step(annealing_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
time_end = time.time()
print("done with annealing in {} minutes\n".format((time_end-time_start)/60.))
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_annealing.pdb'.format(prefix),'w'))
simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
write_full_pdb('{}_post_annealing.pdb'.format(prefix),simulation.topology,positions)
simulation.reporters.pop()


# === Equilibration ===
print('\n=== Equilibrating ===')
simulation.reporters.append(app.statedatareporter.StateDataReporter('{}_thermo_equilibration.out'.format(prefix), thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
print('\nupdating temperature to {} '.format(temperature))
simulation.integrator.setTemperature(temperature) #assumes Langevin
if run_npt:
    simulation.context.setParameter(mm.MonteCarloBarostat.Temperature(), temperature)
time_start = time.time()
simulation.step(equilibration_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
time_end = time.time()
print("done with equilibration in {} minutes\n".format((time_end-time_start)/60.))
simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
write_full_pdb('{}_post_equilibration.pdb'.format(prefix),simulation.topology,positions)
#app.pdbfile.PDBFile.writeHeader(simulation.topology,open('{}_post_equilibration.pdb'.format(prefix),'w'))
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_equilibration.pdb'.format(prefix),'a'))
#app.pdbfile.PDBFile.writeFooter(simulation.topology,open('{}_post_equilibration.pdb'.format(prefix),'a'))

# === Run and save output, with connectivity ===
print('\n=== Production run ===')
simulation.reporters.append(app.statedatareporter.StateDataReporter('{}_thermo_production.out'.format(prefix), thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
time_start = time.time()
if run_npt:
    simulation.reporters.append(app.dcdreporter.DCDReporter('{}_output.dcd'.format(prefix), dcd_report_freq))
else:#nvt
    simulation.reporters.append(app.dcdreporter.DCDReporter('{}_output.dcd'.format(prefix), dcd_report_freq, enforcePeriodicBox=False))

while production_steps > 0:
    nsteps = min( production_steps, 10*dcd_report_freq )
    simulation.step(nsteps)
    production_steps -= nsteps
    # save checkpoints
    simulation.saveCheckpoint('{}_checkpoint.chk'.format(prefix))
    simulation.saveState('{}_checkpoint.xml'.format(prefix))

#simulation.step(production_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
box = simulation.context.getState().getPeriodicBoxVectors(asNumpy=True)
time_end = time.time()
print("done with production in {} minutes\n".format((time_end-time_start)/60.))


simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
write_full_pdb('{}_post_production.pdb'.format(prefix),simulation.topology,positions)
#app.pdbfile.PDBFile.writeHeader(simulation.topology,open('{}_post_production.pdb'.format(prefix),'w'))
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_production.pdb'.format(prefix),'a'))
#app.pdbfile.PDBFile.writeFooter(simulation.topology,open('{}_post_production.pdb'.format(prefix),'a'))




