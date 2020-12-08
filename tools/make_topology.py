# Script to read in single-molecule structures, prepare pdb box, and relevant topology output
# Two operation modes:
#   1) write out a full pdb with connectivity specific to this system
#   2) write out an easily-modifiable .txt file that this function can then use to generate an openmm topology object
#
# example usage:
#   python make_topology.py -L 10.0 -pkmldir /home/kshen/openmm/Structures/packmol -structlib ../../structlib/ -m polystyrene/PS10_a1.pdb 10 -m polystyrene/PS10_a2.pdb 10 --fflist ../TrappeUA_Styrene_Gromos.xml

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

def write_full_pdb(filename,topology,positions):
    app.pdbfile.PDBFile.writeHeader(topology,open(filename,'w'))
    app.pdbfile.PDBFile.writeModel(topology,positions,open(filename,'a'))
    app.pdbfile.PDBFile.writeFooter(topology,open(filename,'a'))

# === Simulation Options ===
prefix = 'packing'
pressure        = 1.0  # bar
temperature     = 773. # kelvin
barostatfreq    = 25
#box_L           = 10.  #nm
run_npt         = True
use_gpu         = True

#nonbonded_method = app.LJPME #set below
nonbonded_cutoff = 1.*unit.nanometer
ewald_tol = 0.0001
friction = 1./unit.picosecond
dt = 0.002 #unit.picosecond

dcd_report_freq = int(1000./dt)
thermo_report_freq = int(1./dt)
equilibration_steps = int(100./dt)
production_steps = int(1.e2/dt)

_C1 = app.element.Element(201, "C1", "CC", 13.018*unit.dalton)
_C2 = app.element.Element(202, "C2", "CT", 13.018*unit.dalton)


# === Parse ===
parser = ap.ArgumentParser(description="Generate a simulation's topology file from single-molecule pdbs")
parser.add_argument('-pkmldir', type=str, default='.', help = "directory to packmol executable")
parser.add_argument('-structlib', type=str, default='.', help = "directory to library of structures")
parser.add_argument('-L', type=float, default=10., help = '(cubic) box L in nm')
parser.add_argument('-m','--molpair', action='append', nargs=2, help = "pairs: molecule_file_name molecule_number")
parser.add_argument('-top','--topfile', type=str, default='', help = "optionally, enter a text file that has the molecule file name and molecule #'s specified")
parser.add_argument('-ff','--fflist', type=str, nargs='+', help = "string of path to ff files")
parser.add_argument('-PME', action='store_true', help = 'LJPME is default, toggle to enable PME')
parser.add_argument('-tail', action='store_true', help = 'tail correction is false by default, toggle to turn on')
args = parser.parse_args()

box_L = args.L

if args.topfile is not '':
    raise ValueError('Not implmeneted yet')
    with open(args.topfile,'r') as f:
        sys_descrip = []
        for line in f:
            tmp = line.split()
            sys_descrip.append( [(tmp[0], int(tmp[1]))] )
else:
    with open('topology.txt','w') as f:
        for entry in args.molpair:
            f.write('{}\t{}\n'.format(entry[0],entry[1]))
    sys_descrip = [ (entry[0], int(entry[1])) for entry in args.molpair ]

sys = []
nMols = []
chainPDBs = []
for entry in sys_descrip:
    pdb = app.PDBFile('{}{}'.format(args.structlib,entry[0]))
    sys.append( (pdb.topology, entry[1]) )
    chainPDBs.append(os.path.join(args.structlib,entry[0]))
    nMols.append(entry[1])

#ff_list = ['TrappeUA_Styrene_Gromos.xml','tip4pew.xml']
ff_list = args.fflist


if args.PME:
    print('using PME with tail correction {}'.format(args.tail))
    nonbonded_method = app.PME
    use_tail = args.tail
else: 
    print('using LJPME, tail correction automatically off')
    nonbonded_method = app.LJPME
    use_tail = False

# === Make topology ===
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

#top_xml = mm.openmm.XmlSerializer.serialize(top)
#with open('proposed_top.xml','w') as f:
#    f.write(top_xml)


# === Packmol ===
def PackMol(nMols, chainPDBs, x, y, z, pdbMix = '{}_packmol.pdb'.format(prefix), logFile = '{}_packmol.log'.format(prefix),mixFile = '{}_packmol.inp'.format(prefix)): 
    """nMols: list of number of molecules for each species
    chainPDBs: list of pdbs of each species
    x,y,z: box dimensions in nm"""
    #convert to angstrom and subtrac 1 angstrom
    x,y,z = (x*10., y*10., z*10.)
    s = """tolerance 2.0
filetype pdb
output {}\n""".format(pdbMix)

    for i, nMol in enumerate(nMols):
        if nMol > 0:
            s += """structure {pdb}
\tnumber {n}
\tresnumbers 2
\tinside box 1 1 1 {x} {y} {z}
\tend structure\n""".format(n = nMol, pdb=chainPDBs[i], x=x-1, y=y-1, z=z-1)
    file = open(mixFile,'w')
    file.write(s)
    file.close()

    print('Packing molecules ...')
    os.system('{}/packmol < {} > {}'.format(args.pkmldir, mixFile, logFile))
    finished = False
    while not finished:
        f = open(logFile,'r')
        if 'Success' in f.read():
            finished = True
        time.sleep(1) 
    return pdbMix

# === Set up short simulation to equilibrate/pack ===
forcefield = app.ForceField(*ff_list)
unmatched_residues = forcefield.getUnmatchedResidues(top)
print('\n=== Unmatched residues ===\n')
print(unmatched_residues)


print('\n=== Periodic Box ===')
periodic_box_vectors = [[box_L,0.,0.],[0.,box_L,0.],[0.,0.,box_L]]
print(periodic_box_vectors)
top.setPeriodicBoxVectors(periodic_box_vectors*unit.nanometer)


print('\n=== Making System ===')
system = forcefield.createSystem(top, 
                                nonbondedMethod = nonbonded_method,
                                nonbondedCutoff = nonbonded_cutoff, 
                                ewaldErrorTolerance=ewald_tol, 
                                rigidWater=True, 
                                constraints=None)


barostat = mm.MonteCarloBarostat( pressure*unit.bar, temperature*unit.kelvin, barostatfreq )
#barostat = mm.MonteCarloBarostat( pressure, temperature, barostatfreq )
if run_npt:
    system.addForce(barostat)


integrator = mm.LangevinIntegrator(temperature*unit.kelvin, friction, dt*unit.picosecond)
if use_gpu:
    platform = mm.Platform.getPlatformByName('OpenCL')
    properties = {'DeviceIndex':'0', 'Precision':'mixed'}
else:
    platform = mm.Platform.getPlatformByName('CPU')
    properties = {'Threads': '1'}

simulation = app.Simulation( top, system, integrator, platform, properties )


# === Get initial configuration === 
out_pdb_filename = PackMol(nMols, chainPDBs, box_L, box_L, box_L) 
sys_pdb = app.PDBFile(out_pdb_filename)
positions = sys_pdb.positions


simulation.context.setPositions(positions)
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_packmol.pdb'.format(prefix),'w'))
write_full_pdb('{}_packmol.pdb'.format(prefix),simulation.topology,positions)
simulation.context.setPeriodicBoxVectors(periodic_box_vectors[0],periodic_box_vectors[1],periodic_box_vectors[2]) # Set the periodic box vectors


simulation.context.applyConstraints(1e-8)


# By default PME turns on tail correction. Manually turn off if requested
print('setting tail correction to {}'.format(use_tail))
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
simulation.minimizeEnergy(tolerance=0.01*unit.kilojoules_per_mole,maxIterations=1000000)
print('post minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
time_end = time.time()
print("done with minimization in {} minutes\n".format((time_end-time_start)/60.))
positions = simulation.context.getState(getPositions=True).getPositions()       
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_minimization.pdb'.format(prefix),'w'))
simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
write_full_pdb('{}_post_minimization.pdb'.format(prefix),simulation.topology,positions)

# === Setup Run ===
simulation.reporters.append(app.statedatareporter.StateDataReporter('{}_thermo.out'.format(prefix), thermo_report_freq, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True))

# === Equilibration ===
time_start = time.time()
simulation.step(equilibration_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
time_end = time.time()
print("done with equilibration in {} minutes\n".format((time_end-time_start)/60.))
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_post_equilibration.pdb'.format(prefix),'w'))
simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
write_full_pdb('{}_post_equilibration.pdb'.format(prefix),simulation.topology,positions)

# === Run and save output, with connectivity ===
time_start = time.time()
simulation.reporters.append(app.dcdreporter.DCDReporter('{}_output.dcd'.format(prefix), dcd_report_freq))
simulation.step(production_steps)
positions = simulation.context.getState(getPositions=True).getPositions()       
box = simulation.context.getState().getPeriodicBoxVectors(asNumpy=True)
time_end = time.time()
print("done with production in {} minutes\n".format((time_end-time_start)/60.))


simulation.topology.setPeriodicBoxVectors( simulation.context.getState().getPeriodicBoxVectors() )
#app.pdbfile.PDBFile.writeHeader(simulation.topology,open('{}_proposal.pdb'.format(prefix),'w'))
#app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('{}_proposal.pdb'.format(prefix),'a'))
#app.pdbfile.PDBFile.writeFooter(simulation.topology,open('{}_proposal.pdb'.format(prefix),'a'))
write_full_pdb('{}_proposal.pdb'.format(prefix),simulation.topology,positions)
#np.savetxt( 'proposed_box.txt', np.vstack([np.diag(box),box] ) )
np.savetxt( '{}_box.txt'.format(prefix), box )




