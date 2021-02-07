#Replicate a box
# Version 1: January 2021, Dan Sun and Kevin Shen
# Right now, assumes replicating into a cube
#
# example:
# python -i replicate.py run0_output.dcd 8 -padding 2 -topfile packing_packmol.pdb
# python -i replicate.py packing_packmol.pdb 8 -padding 2
#


import numpy as np
import mdtraj
import argparse as ap

parser = ap.ArgumentParser(description='Replicate a cell')
parser.add_argument('coordfile',type=str,help='coordinate file name')
parser.add_argument('nimages',type=int,help='number of total unit cells to have at the end')
parser.add_argument('-topfile',default=None,type=str,help='topology file, if needed (e.g. for dcd)')
parser.add_argument('-padding',type=int,default=0, help='put additional padding between cells')
parser.add_argument('-L',type=float,default=None, help='unit cell size, if not included in coordinate file')
args = parser.parse_args()

suffix = args.coordfile.split('.')[-1]
prefix = '.'.join(args.coordfile.split('.')[:-1])
dim = 3


# initializing
print('... Initializing ...')
if args.topfile is not None:
    traj = mdtraj.load(args.coordfile, top=args.topfile)
else:
    traj = mdtraj.load(args.coordfile)

xyz0 = traj.xyz[0]
if args.L is None:
    boxL = traj.unitcell_lengths[0]
else:
    boxL = [args.L, args.L, args.L]
boxL = np.array(boxL)
padded_box = boxL + np.array(3*[args.padding])

# replicate the topology
print('... Replicating topology ...')
new_top = traj.top.copy()
for ii in np.arange(1,args.nimages):
    print(' adding topology for image {}'.format(ii))
    new_top = new_top.join( traj.top.copy() )
    

# replicate the coordinates
# assuming that final cell is also cubic, for the time being
# choose random cells to fill
# in the future, can also intelligently figure out non-cubic cells, or allow user to input the final tiled shape
print('... Replicating coordinates ...')
cell_scaling = 1
cell_scaling = np.int( np.ceil( args.nimages**(1.0/dim) ) )
print('cell scaling: {}'.format(cell_scaling))

indices = np.random.choice( np.arange(cell_scaling**3), size=args.nimages, replace=False )
xyz = np.zeros([traj.n_atoms * args.nimages, 3])

xyz_rand = traj.xyz[ np.random.choice(np.arange(traj.n_frames), size=args.nimages), :,: ] #choose random frames, with replacement (i.e. if fewer frames in trajectry than desired images)
for ii,index in enumerate(indices):
    index_tuple = np.unravel_index(index,3*[cell_scaling])
    print('allocating coordinates to image {} with cell index {}'.format(index,index_tuple))

    xyz_tmp = xyz_rand[ii,:,:]
    xyz[ii*traj.n_atoms:(ii+1)*traj.n_atoms,:] = xyz_tmp + padded_box * index_tuple 


# save
print('... Saving ...')
new_traj = mdtraj.Trajectory(np.array([xyz]),new_top,unitcell_lengths=padded_box,unitcell_angles=[90.,90.,90.])
new_traj.save(prefix+'_{}image.{}'.format(args.nimages,suffix))


