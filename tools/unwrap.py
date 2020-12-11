import argparse
import mdtraj as md
import numpy as np


def unwrap_traj(traj, top, method='npt'):

    # import wrapped trajectory
    traj_w = md.load(traj, top=top)
    xyz_w = traj_w.xyz

    # get prefix from traj filename
    prefix = '.'.join(traj.split('.')[:-1])

    # initialize unwrapped positions
    xyz_uw = np.empty([traj_w.n_frames, traj_w.topology.n_atoms, 3])

    # frame 0 in unwrapped trajectory is same as the wrapped trajectory's
    xyz_uw[0] = xyz_w[0]

    # unwrap trajectories
    if method == 'nvt':
        for i in range(1, len(xyz_uw)):
            uc_vecs = traj_w.unitcell_vectors[i]
            uc_vecs_inv = np.linalg.inv(uc_vecs)
            xyz_uw[i] = xyz_w[i] - np.floor((xyz_w[i] - xyz_uw[i-1]).dot(uc_vecs_inv) + np.array([0.5, 0.5, 0.5])).dot(uc_vecs)
    elif method == 'npt':
        for i in range(1, len(xyz_uw)):
            uc_vecs = traj_w.unitcell_vectors[i]
            uc_vecs_inv = np.linalg.inv(uc_vecs)
            xyz_uw[i] = xyz_uw[i-1] + (xyz_w[i] - xyz_w[i-1]) - np.floor((xyz_w[i] - xyz_w[i-1]).dot(uc_vecs_inv) + np.array([0.5, 0.5, 0.5])).dot(uc_vecs)
    else:
        raise ValueError("method must by 'nvt' or 'npt'")

    # create new unwrapped trajectory
    traj_uw = md.Trajectory(xyz_uw, traj_w.topology,
                            unitcell_lengths=traj_w.unitcell_lengths, unitcell_angles=traj_w.unitcell_angles)

    # save unwrapped trajectory to file
    traj_uw.save("{}_uw.dcd".format(prefix))

    # return unwrapped trajectory
    return traj_uw


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('traj', type=str, help="trajectory file")
    parser.add_argument('top', type=str, help="topology file")
    parser.add_argument('--method', default='nvt', type=str, help="method of unwrapping")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    unwrap_traj(args.traj, args.top, method=args.method)
