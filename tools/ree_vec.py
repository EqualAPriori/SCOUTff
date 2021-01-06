"""
Wishlist:
- sliding window averaging for better statistics
- choosing what window to fit over
"""
import mdtraj
import numpy as np
from scipy.optimize import curve_fit
import argparse as ap

def parse_arguments():
    parser = ap.ArgumentParser(description="calculate msd statistics")
    parser.add_argument('traj', type=str, help = "trajectory file")
    parser.add_argument('top', type=str, help = "topology file")
    parser.add_argument('-stride', type=int, default=1, help = "stride")
    parser.add_argument('-dt', type=float, default=0.1, help = "dt (between frames, before striding)")
    parser.add_argument('-frac', type=float, default=0.1, help = "fraction of data to do Ree_vec correlation fit")
    args = parser.parse_args()
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    t  = mdtraj.load(args.traj, top=args.top, stride=args.stride)
    dt = args.dt * args.stride
    frac = args.frac


    indices = np.zeros([t.topology.n_chains,2],dtype='int')

    for ic,c in enumerate(t.topology.chains):
        atom_list = list(c.atoms)
        indices[ic,:] = [atom_list[0].index, atom_list[-1].index]


    Ree_vec = t.xyz[:,indices[:,1],:] - t.xyz[:,indices[:,0],:]
    Ree_vec0 = Ree_vec[0,:,:]

    Ree_dot = np.sum(Ree_vec * Ree_vec0[None,:,:],2)
    Ree_dot_avg = np.mean(Ree_dot,1)

    #Fit exponential only to the first 10% of data
    def my_exp(t, a,tau):
        return a*np.exp(-t/tau)

    fracs = [0.01,0.025,0.05,0.1,0.2,0.3,0.5]
    Cs = np.zeros(len(fracs))
    taus = np.zeros(len(fracs))

    for ifrac,frac in enumerate(fracs):
        print('...using first {} of data'.format(frac))
        time = np.arange(t.n_frames)*dt
        ind_max = int(np.floor(t.n_frames*frac))

        popt, pcov = curve_fit(my_exp, time[:ind_max], Ree_dot_avg[:ind_max])

        print('{} * exp(-t/{})'.format(popt[0],popt[1]))
        print('absolute_sigma=False pcov:\n{}'.format(pcov))
    
        Cs[ifrac] = popt[0]
        taus[ifrac] = popt[1]

    data = np.vstack([fracs,Cs,taus]).T
    print(data)
    np.savetxt('ree_rvec_corr.dat',data,header='fraction_of_data\tC\ttau;\t C*exp(-t/tau)')

