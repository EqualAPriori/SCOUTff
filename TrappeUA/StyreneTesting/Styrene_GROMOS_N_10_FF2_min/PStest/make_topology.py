# Script to read in single-molecule structures, prepare box, and relevant topology output
# Two operation modes:
#   1) write out a full pdb with connectivity specific to this system
#   2) write out an easily-modifiable .txt file that this function can then use to generate an openmm topology object
#
# example usage:
#   python make_topology.py -lib ../../../../structlib/ -m C12.pdb 10 -m cl.pdb 5 -m polystyrene/PS10_a1.pdb 10
#   python make_topology.py -lib ../../../../structlib/ -m polystyrene/PS10_a1.pdb 10 -m polystyrene/PS10_a2.pdb 10 --fflist ../TrappeUA_Styrene_Gromos.xml ../tip4pew.xml
#   python make_topology.py -lib ../../../../structlib/ -m polystyrene/PS10_a1.pdb 10 -m polystyrene/PS10_a2.pdb 10 --fflist ../TrappeUA_Styrene_Gromos.xml
#
# key variables:
#   sys_descrip: [ ('mol1filename',num_mol1), ('mol2filename',num_mol2), ... ]
#   sys: [ (mol1_ommtopology,num_mol1), (mol2_ommtopology,num_mol2), ... ]
#   top: the openmm topology with desired # molecules
#
# Todo:
#   1) add packmol engine
#   2) add arguments to specify library for xml forcefields... maybe can also specify with a text file, kind of like specification of topology?
#   3) barostat giving errors right now...
#

from simtk import unit
from simtk.openmm import app
import simtk.openmm as mm
import numpy as np
import argparse as ap
import time
import os
from subprocess import call
#def molpair(arg):
#    # For simplity, assume arg is a pair of integers separated by a comma. 
#    return [arg[0], int(arg[1])]

parser = ap.ArgumentParser(description="Generate a simulation's topology file from single-molecule pdbs")
parser.add_argument('-lib', type=str, default='.', help = "directory to library of structures")
parser.add_argument('-m','--molpair', action='append', nargs=2, help = "pairs: molecule_file_name molecule_number")
parser.add_argument('-top','--topfile', type=str, default='', help = "optionally, enter a text file that has the molecule file name and molecule #'s specified")
parser.add_argument('-ff','--fflist', type=str, nargs='+', help = "string of path to ff files")
args = parser.parse_args()

# === Parse ===
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
    pdb = app.PDBFile('{}{}'.format(args.lib,entry[0]))
    sys.append( (pdb.topology, entry[1]) )
    chainPDBs.append(os.path.join(args.lib,entry[0]))
    nMols.append(entry[1])

#ff_list = ['TrappeUA_Styrene_Gromos.xml','tip4pew.xml']
ff_list = args.fflist

# === Make topology ===
top = app.topology.Topology()
for ii,mol in enumerate(sys):
    moltop, nummol = mol
    for jj in range(nummol):
        atoms_in_top = []
        for c in moltop.chains():
            chain = top.addChain()
            for r in c.residues():
                residue = top.addResidue("custom",chain)
                for a in enumerate(r.atoms()):
                    atom = top.addAtom(a[1].name, a[1].element, residue)
                    atoms_in_top.append(atom)
                    print(atom)
        for bond in moltop.bonds():
            bi1,bi2 = bond[0].index, bond[1].index
            top.addBond( atoms_in_top[bi1], atoms_in_top[bi2] ) 

print(top)

#test = mm.openmm.XmlSerializer.serialize(top)
#with open('top.xml','w') as f:
#    f.write(test)


# === Packmol ===
def PackMol(nMols, chainPDBs, x, y, z, pdbMix = 'initial.pdb', logFile = 'packmol.log',mixFile = 'mix.inp'): 
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
    os.system('packmol < {} > {}'.format(mixFile, logFile))
    finished = False
    while not finished:
        f = open(logFile,'r')
        if 'Success' in f.read():
            finished = True
        time.sleep(1) 
    return pdbMix

# === Set up short simulation to equilibrate/pack ===
pressure        = 1.0  # bar
temperature     = 473. # kelvin
barostatfreq    = 25.
box_L            = 10.  #nm
run_npt         = True
use_gpu         = True



forcefield = app.ForceField(*ff_list)
unmatched_residues = forcefield.getUnmatchedResidues(top)
print('Unmatched\n')
print(unmatched_residues)


periodic_box_vectors = [[box_L,0.,0.],[0.,box_L,0.],[0.,0.,box_L]]
top.setPeriodicBoxVectors(periodic_box_vectors*unit.nanometer)


system = forcefield.createSystem(top, 
                                nonbondedMethod=app.LJPME,
                                nonbondedCutoff=1.*unit.nanometer, 
                                ewaldErrorTolerance=0.0001, 
                                rigidWater=True, 
                                constraints=None)


integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1./unit.picosecond, 0.002*unit.picoseconds)
#barostat = mm.MonteCarloBarostat( pressure*unit.bar, temperature*unit.kelvin, barostatfreq )
#barostat = mm.MonteCarloBarostat( pressure, temperature, barostatfreq )
#if run_npt:
#    system.addForce(barostat)

if use_gpu:
    platform = mm.Platform.getPlatformByName('OpenCL')
    properties = {'DeviceIndex':'1', 'Precision':'mixed'}
else:
    platform = Platform.getPlatformByName('CPU')
    properties = {'Threads': '1'}

simulation = app.Simulation( top, system, integrator, platform, properties )

#UPDATE to use packmol positions
out_pdb_filename = PackMol(nMols, chainPDBs, box_L, box_L, box_L) 
sys_pdb = app.PDBFile(out_pdb_filename)
positions = sys_pdb.positions
#positions = np.random.random((top.getNumAtoms(),3))*box_L

simulation.context.setPositions(positions)
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('initial.pdb','w'))
simulation.context.setPeriodicBoxVectors(periodic_box_vectors[0],periodic_box_vectors[1],periodic_box_vectors[2]) # Set the periodic box vectors


simulation.context.applyConstraints(1e-8)

time_start = time.time()
print('pre minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy(tolerance=0.01*unit.kilojoules_per_mole,maxIterations=1000000)
print('post minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
time_end = time.time()
print("done with minimization in {} mininutes\n".format((time_end-time_start)/60.))
simulation.reporters.append(app.dcdreporter.DCDReporter('npt.dcd', 1000))
simulation.reporters.append(app.statedatareporter.StateDataReporter('thermo.out', 1000, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True))


# === Run and save output, with connectivity ===
simulation.step(5000)
positions = simulation.context.getState(getPositions=True).getPositions()       
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('short_run.pdb','w'))

simulation.step(5E6)
positions = simulation.context.getState(getPositions=True).getPositions()       
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('longer_run.pdb','w'))




