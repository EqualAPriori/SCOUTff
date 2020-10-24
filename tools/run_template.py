# Script to read in single-molecule structures, prepare pdb box, and relevant topology output
# Two operation modes:
#   1) write out a full pdb with connectivity specific to this system
#   2) write out an easily-modifiable .txt file that this function can then use to generate an openmm topology object
#
# example usage:
#   python run.py -L 10.0 --fflist ./TrappeUA_Styrene_Gromos.xml
# 
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

# === Simulation Options ===
pressure        = 1.0  # bar
temperature     = __TEMP__ # kelvin
temperature_anneal = __TEMPANNEAL__ #kelvin
barostatfreq    = 25
run_npt         = True
use_gpu         = True

nonbonded_method = app.LJPME
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
parser.add_argument('-ff','--fflist', type=str, nargs='+', help = "string of path to ff files")
args = parser.parse_args()

prefix = args.prefix
if args.L is not None:
    print('CAUTION: (cubic) box L {} entered, overriding boxfile {}'.format(args.L, args.boxfile))
    box_L = args.L
    periodic_box_vectors = [[box_L,0.,0.],[0.,box_L,0.],[0.,0.,box_L]]
else:
    periodic_box_vectors = np.loadtxt( args.boxfile ) 

ff_list = args.fflist
pdb = app.PDBFile 

sys_pdb = app.PDBFile(args.initpdb)
top = sys_pdb.topology

# === Set up simulation ===
forcefield = app.ForceField(*ff_list)
unmatched_residues = forcefield.getUnmatchedResidues(top)
print('\n=== Unmatched residues ===\n')
print(unmatched_residues)


print('\n=== Periodic Box ===')
print(periodic_box_vectors)
top.setPeriodicBoxVectors(periodic_box_vectors*unit.nanometer)


print('\n=== Making System ===')
system = forcefield.createSystem(sys_pdb.topology, 
                                nonbondedMethod = nonbonded_method,
                                nonbondedCutoff = nonbonded_cutoff, 
                                ewaldErrorTolerance=ewald_tol, 
                                rigidWater=True, 
                                constraints=None)


integrator = mm.LangevinIntegrator(temperature_anneal*unit.kelvin, friction, dt*unit.picosecond)
barostat = mm.MonteCarloBarostat( pressure*unit.bar, temperature_anneal*unit.kelvin, barostatfreq )
#barostat = mm.MonteCarloBarostat( pressure, temperature, barostatfreq )
if run_npt:
    system.addForce(barostat)

if use_gpu:
    platform = mm.Platform.getPlatformByName('OpenCL')
    properties = {'DeviceIndex':'1', 'Precision':'mixed'}
else:
    platform = mm.Platform.getPlatformByName('CPU')
    properties = {'Threads': '1'}

simulation = app.Simulation( top, system, integrator, platform, properties )


# === Get initial configuration === 
positions = sys_pdb.positions


simulation.context.setPositions(positions)
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_initial.pdb'.format(prefix),'w'))
simulation.context.setPeriodicBoxVectors(periodic_box_vectors[0],periodic_box_vectors[1],periodic_box_vectors[2]) # Set the periodic box vectors


simulation.context.applyConstraints(1e-8)


# === Minimize Energy ===
print('\n=== Minimizing ===')
time_start = time.time()
print('pre minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy(tolerance=0.01*unit.kilojoules_per_mole,maxIterations=1000000)
print('post minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
time_end = time.time()
print("done with minimization in {} minutes\n".format((time_end-time_start)/60.))
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_minimization.pdb'.format(prefix),'w'))


# === Setup Runs ===

# === Annealing ===
print('\n=== Annealing ===')
simulation.reporters.append(app.statedatareporter.StateDataReporter('{}_thermo_annealing.out'.format(prefix), thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
print('\nannealing at temperature to {} '.format(temperature_anneal))
time_start = time.time()
simulation.step(equilibration_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
time_end = time.time()
print("done with annealing in {} minutes\n".format((time_end-time_start)/60.))
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_annealing.pdb'.format(prefix),'w'))
simulation.reporters.pop()


# === Equilibration ===
print('\n=== Equilibrating ===')
simulation.reporters.append(app.statedatareporter.StateDataReporter('{}_thermo_equilibration.out'.format(prefix), thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
print('\nupdating temperature to {} '.format(temperature))
simulation.integrator.setTemperature(temperature) #assumes Langevin
simulation.context.setParameter(mm.MonteCarloBarostat.Temperature(), temperature)
time_start = time.time()
simulation.step(equilibration_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
time_end = time.time()
print("done with equilibration in {} minutes\n".format((time_end-time_start)/60.))
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_equilibration.pdb'.format(prefix),'w'))

# === Run and save output, with connectivity ===
print('\n=== Production run ===')
simulation.reporters.append(app.statedatareporter.StateDataReporter('{}_thermo_production.out'.format(prefix), thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
time_start = time.time()
simulation.reporters.append(app.dcdreporter.DCDReporter('{}_output.dcd'.format(prefix), dcd_report_freq))
simulation.step(production_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
box = simulation.context.getState().getPeriodicBoxVectors(asNumpy=True)
time_end = time.time()
print("done with production in {} minutes\n".format((time_end-time_start)/60.))


simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
app.pdbfile.PDBFile.writeHeader(simulation.topology,open('{}_post_production.pdb'.format(prefix),'w'))
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_production.pdb'.format(prefix),'a'))
app.pdbfile.PDBFile.writeFooter(simulation.topology,open('{}_post_production.pdb'.format(prefix),'a'))




