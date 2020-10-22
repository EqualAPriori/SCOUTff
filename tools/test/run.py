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
temperature     = 473. # kelvin
barostatfreq    = 25
#box_L           = 10.  #nm
run_npt         = True
use_gpu         = True

nonbonded_method = app.LJPME
nonbonded_cutoff = 1.*unit.nanometer
ewald_tol = 0.0001
friction = 1./unit.picosecond
dt = 0.002 #unit.picosecond

dcd_report_freq = int(1000./dt)
thermo_report_freq = int(1./dt)

#Annealing
equilibration_steps = int(5.0e2/dt)
production_steps = int(1.0e2/dt)
_C1 = app.element.Element(201, "C1", "CC", 13.018*unit.dalton)
_C2 = app.element.Element(202, "C2", "CT", 13.018*unit.dalton)


# === Parse ===
parser = ap.ArgumentParser(description="Generate a simulation's topology file from single-molecule pdbs")
parser.add_argument('-initpdb', type=str, default='proposed_initial_packing.pdb', help = "initial config. MUST have connectivity.")
parser.add_argument('-L', type=float, default=10., help = '(cubic) box L in nm')
parser.add_argument('-ff','--fflist', type=str, nargs='+', help = "string of path to ff files")
args = parser.parse_args()

box_L = args.L
ff_list = args.fflist
pdb = app.PDBFile 

sys_pdb = app.PDBFile(args.initpdb)
top = sys_pdb.topology

# === Set up simulation ===
forcefield = app.ForceField(*ff_list)
unmatched_residues = forcefield.getUnmatchedResidues(top)
print('Unmatched\n')
print(unmatched_residues)


periodic_box_vectors = [[box_L,0.,0.],[0.,box_L,0.],[0.,0.,box_L]]
top.setPeriodicBoxVectors(periodic_box_vectors*unit.nanometer)


system = forcefield.createSystem(sys_pdb.topology, 
                                nonbondedMethod = nonbonded_method,
                                nonbondedCutoff = nonbonded_cutoff, 
                                ewaldErrorTolerance=ewald_tol, 
                                rigidWater=True, 
                                constraints=None)


integrator = mm.LangevinIntegrator(temperature*unit.kelvin, friction, dt*unit.picosecond)
barostat = mm.MonteCarloBarostat( pressure*unit.bar, temperature*unit.kelvin, barostatfreq )
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
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('initial.pdb','w'))
simulation.context.setPeriodicBoxVectors(periodic_box_vectors[0],periodic_box_vectors[1],periodic_box_vectors[2]) # Set the periodic box vectors


simulation.context.applyConstraints(1e-8)


# === Minimize Energy ===
time_start = time.time()
print('pre minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy(tolerance=0.01*unit.kilojoules_per_mole,maxIterations=1000000)
print('post minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
time_end = time.time()
print("done with minimization in {} minutes\n".format((time_end-time_start)/60.))
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('post_minimization.pdb','w'))


# === Setup Run ===
simulation.reporters.append(app.statedatareporter.StateDataReporter('thermo.out', thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True))

# === Equilibration ===
time_start = time.time()
simulation.step(equilibration_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
time_end = time.time()
print("done with equilibration in {} minutes\n".format((time_end-time_start)/60.))
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('post_equilibration.pdb','w'))

# === Run and save output, with connectivity ===
time_start = time.time()
simulation.reporters.append(app.dcdreporter.DCDReporter('output.dcd', dcd_report_freq))
simulation.step(production_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
time_end = time.time()
print("done with production in {} minutes\n".format((time_end-time_start)/60.))

app.pdbfile.PDBFile.writeHeader(simulation.topology,open('post_production.pdb','w'))
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('post_production.pdb','a'))
app.pdbfile.PDBFile.writeFooter(simulation.topology,open('post_production.pdb','a'))




