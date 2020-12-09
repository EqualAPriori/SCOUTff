#Wishlist:
#1) read in a file with mass definitions for center of mass calculation
#2) same, but for chain names (since pdb doesn't define chain name)
#
#Example usage: 
#   python ~/mylib/SCOUTff/tools/com.py run3_output.dcd run3_initial.pdb
#   python ~/mylib/SCOUTff/tools/com.py run3_output.dcd run3_initial.pdb -com
import mdtraj
import numpy as np
import argparse as ap

parser = ap.ArgumentParser(description="calculate centroid or center of mass")
parser.add_argument('traj', type=str, help = "trajectory file")
parser.add_argument('top', type=str, help = "topology file")
parser.add_argument('-com', action='store_true', help = "toggle on center of mass calculation, using elemental masses")
args = parser.parse_args()

prefix = '.'.join(args.traj.split('.')[:-1])
label = 'com' if args.com else 'centroid'
traj = mdtraj.load(args.traj, top=args.top)

#Create new topology
bead_dict = {}
num_bead_types = 0
chain_name = 'my_chain' #chain_name is not stored anywhere
new_topology = mdtraj.Topology()
for chain in traj.topology.chains:
    masses_in_chain = np.array([atom.element.mass for residue in chain.residues for atom in residue.atoms])
    mass_chain = masses_in_chain.sum()

    new_chain = new_topology.add_chain()
    
    new_res = new_topology.add_residue(chain_name,new_chain)
    if chain_name not in bead_dict:
        bead_dict[chain_name] = mdtraj.element.Element(200+num_bead_types,chain_name,'A'+str(num_bead_types),mass_chain,1.0)
        num_bead_types += 1
    new_bead = new_topology.add_atom(chain_name,bead_dict[chain_name],new_res)
   
#Get trajectory
new_xyz = np.zeros([traj.n_frames,new_topology.n_chains,3])
if args.com:
    print('Currently use inferred bead masses to calculate com')

for ic,chain in enumerate(traj.topology.chains):
    atom_indices = [atom.index for residue in chain.residues for atom in residue.atoms]
    if args.com:
        masses_in_chain = np.array([atom.element.mass for residue in chain.residues for atom in residue.atoms])
    else:
        masses_in_chain = np.ones(chain.n_atoms)

    traj_slice = traj.atom_slice(atom_indices)
    xyzs = (traj_slice.xyz*masses_in_chain[None,:,None]).mean(1)
    #xyzs = np.expand_dims(xyzs,1)
    new_xyz[:,ic,:] = xyzs

#Save out
new_traj = mdtraj.Trajectory(new_xyz,new_topology,unitcell_lengths = traj.unitcell_lengths, unitcell_angles = traj.unitcell_angles)
new_traj.save('{}_{}.dcd'.format(prefix,label))
new_traj[0].save('{}_{}.pdb'.format(prefix,label))


