"""
=== 2020.12.14 ===
Things to consider
1) Single polymer chain, sliding window msd analysis, check for correlation times (along t0) for the quantity msd(t0;n*dt) = (r(t0 + n*dt) - r(t0))**2.0
2) Calculate variance & std error of msd across polymer chains (note, Shanbhag 2015 calculated variances across *sliding-window-averaged* msd for each chain)
    We improve upon Shanbhag 2015 a little by only averaging over statistically independent windows, instead of every single window.
3) Bootstrap to calculate confidence interval in diffusion constant
4) Curvature detection / finding suitable sub-interval of msd data to use

Usage: prior to calling this script, need to 1) get com trajectories, and 2) unwrap, if needed.
An example of a full stack of commands is as follows:

python ~/mylib/SCOUTff/tools/com.py run0_output.dcd run0_post_production.pdb &
python ~/mylib/SCOUTff/tools/unwrap.py run0_output_centroid.dcd run0_output_centroid.pdb &
python ~/mylib/SCOUTff/msd.py run0_output_centroid_uw.dcd run0_output_centroid.pdb -stride 100

"""
import mdtraj
import numpy as np
import argparse as ap
from pymbar import timeseries

def parse_arguments():
    parser = ap.ArgumentParser(description="calculate msd statistics")
    parser.add_argument('traj', type=str, help = "trajectory file")
    parser.add_argument('top', type=str, help = "topology file")
    parser.add_argument('-stride', type=int, default=1, help = "stride")
    parser.add_argument('-dt', type=float, default=0.1, help = "dt (between frames, before striding)")
    args = parser.parse_args()
    return parser.parse_args()

#=== Utilities ===
def get_stats(data):
    """
    later, can generalize, to use one column for decorrelating and getting reference indices
    """
    [t0, g, Neff] = timeseries.detectEquilibration(data)
    data_equil = data[t0:]
    indices = timeseries.subsampleCorrelatedData(data_equil, g=g)
    sub_data = data_equil[indices]
    
    avg = sub_data.mean()
    std = sub_data.std()
    err = sub_data.std()/np.sqrt( len(indices) )

    return avg,std,err, t0,g,Neff, sub_data

def calc_msd(_t,frame=0):
    """
    Given one trajectory, also calculate variances about the mean.
    """
    xyz_ref = _t.xyz[frame,:,:]
    disp2 = (_t.xyz - xyz_ref[None,:,:])**2.0
    return np.sum(disp2,2) #n_frames x n_atoms

def calc_msd_windowed(_t,use_mdtraj=False,verbose=False):
    """
    Note, if using mdtraj's rmsd calculation, pre-averages over chains
    first index is initial frame
    second index is time/frame separation
    Wish list:
    - 2020.12.14 program in max initial frame and max time separation
    """
    if use_mdtraj:
        msd_data = np.zeros( [_t.n_frames,_t.n_frames] )
        msd_summary = np.zeros( [_t.n_frames,7] )
        for iframe in range(_t.n_frames):
            if verbose and np.mod(iframe,100) == 0:
                print('frame {}'.format(iframe))
            if iframe == 0:
                msd_data[iframe,:] = mdtraj.rmsd(_t,_t,frame=iframe)[iframe:]**2.0
            else:
                msd_data[iframe,:-iframe] = mdtraj.rmsd(_t,_t,frame=iframe)[iframe:]**2.0
    else:
        msd_data = np.zeros( [_t.n_frames,_t.n_frames,_t.n_atoms] )
        for iframe in range(_t.n_frames):
            if verbose and np.mod(iframe,100) == 0:
                print('frame {}'.format(iframe))
            if iframe == 0:
                msd_data[iframe,:,:] = calc_msd(_t,frame=iframe)[iframe:,:]
            else:
                msd_data[iframe,:-iframe,:] = calc_msd(_t,frame=iframe)[iframe:,:]
     
    return msd_data


def calc_sliding_window_avg(msd_block, n_frames = None):
    """
    takes in [n_frame x n_frame] data, assuming in UL triangular form, 
    where first index is the initial frame index (taken to be the reference)
    and second index is the frame separation from the reference frame

    Note, can be used on per-chain msd_block data, *or* ensemble-averaged msd_block data

    Wishlist: 
    - 2020.12.14 program in max initial frame and max time separation
    """
    if n_frames is None:
        n_frames = np.shape(msd_block)[0]

    msd_summary = np.zeros([n_frames,7])
    for dframe in range(n_frames):
        if np.mod(dframe,100) == 0:
            print('sliding window dframe {}'.format(dframe))
        msd_at_dframe = msd_block[:n_frames-dframe, dframe]

        avg,std,err, t0,g,Neff, sub_data = get_stats(msd_at_dframe)

        msd_summary[dframe,1:] = [avg,std,err, t0,g,Neff] 
        
    msd_summary[:,0] = np.arange(0,t.n_frames)*dt

    return msd_summary
    
def calc_D(time,msd,err=None,n_dim=3):
    """
    Currently uses every time point in the msd-fitting.

    Parameters
    ----------
    time: array-like
        time vector
    msd: array-like
        msd data, preferably already after some suitable averaging process
    err:
        if None, uses the naive result in literature of weighting every time point equally. Otherwise, use errors for Weighted Least Squares.
        in literature, people rely on fitting over much shorter time intervals to get around increasing variance at long times.

    Returns
    -------
    array
        [a,D]
    array
        [a_err, D_err]

    Notes
    -----
    The error estimate here may not be the best -- it's calculated on averaged MSD data already. Maybe bootstrapping is better?

    Wishlist:
    - curvature detection; fitting over a sub-interval. This is maybe an outside sub-routine, after which one feeds in only the desired sub-interval into this function

    References
    ----------
    https://www.stat.colostate.edu/regression_book/chapter8.pdf
        for errors, see eq. 8.2.21, also Eq. 8.2.35 and first two lines of p. 579
    can also use scipy's curve_fit, but here we get error estimate on parameter
    """
    if err is None:
        g = np.ones(time.shape)
    else:
        g = err
        g[g==0] = np.min(g[g>0])/100. #need to filter out zero-error points, e.g. t=0!

    W = np.diag(1/g**2.0)
    X = np.ones([time.size,2])
    X[:,1] = time
    Y = np.reshape(msd,[msd.size,1])

    Cw = np.dot(X.T,np.dot(W,X))
    tmp = np.dot(X.T,np.dot(W,Y))
    beta = np.dot(np.linalg.inv(Cw),tmp) 

    wsse = np.dot(W,((Y-np.dot(X,beta))**2.0)).sum() #weighted sum squared error, Eq. 8.2.7
    wmse = wsse/(Y.size - 2) 
    beta_err = wmse*np.diag(np.linalg.inv(Cw))

    beta = np.ravel(beta)
    beta_err = np.ravel(beta_err)
    print('D = {} +/- {}'.format(beta[1]/2/n_dim, beta_err[1]/2/n_dim))

    return np.array([beta[0],beta[1]/2/n_dim]), np.array([beta_err[0],beta_err[1]/2/n_dim])


def msd_sliding_window(_t):
    """
    DEPRECATED
    Note, this was a preliminary implemenation -- the decorrelation time and statistics are calculated over msd-data where each frame has already been averaged over the chains in the box.
    This thus neglects the chain-level decorrelation time of msd data, and misses the variance.

    Conceptually, there is the question of whether sliding window is appropriate -- i.e. people naively treat a trajectory as independent from t+dt. 
    Intuitively, they should be correlated, reducing the # of independent windows from which one can collect meaningful data.
    """
    msd_data = np.zeros( [_t.n_frames,_t.n_frames] )
    msd_summary = np.zeros( [_t.n_frames,7] )
    for iframe in range(_t.n_frames):
        if np.mod(iframe,100) == 0:
            print('frame {}'.format(iframe))
        if iframe == 0:
            msd_data[iframe,:] = mdtraj.rmsd(_t,_t,frame=iframe)[iframe:]**2.0
        else:
            msd_data[iframe,:-iframe] = mdtraj.rmsd(_t,_t,frame=iframe)[iframe:]**2.0

    for dframe in range(t.n_frames):
        if np.mod(dframe,100) == 0:
            print('dframe {}'.format(dframe))
        msd_at_dframe = msd_data[:_t.n_frames-dframe, dframe]

        avg,std,err, t0,g,Neff, sub_data = get_stats(msd_at_dframe)

        #msd_mean[delta_t] =
        #msd_std[delta_t] =
        #msd_error[delta_t] = 

        msd_summary[dframe,1:] = [avg,std,err, t0,g,Neff] 
        
    msd_summary[:,0] = np.arange(0,_t.n_frames)*dt

    np.savetxt('msd_summary.txt',msd_summary,header="avg\tstd\terr\tt0\tg\tNeff")


#=== Automating MSD calculation ===
if __name__ == '__main__':
    args = parse_arguments()
    t  = mdtraj.load(args.traj, top=args.top, stride=args.stride)
    dt = args.dt * args.stride
    n_frames = t.n_frames
    n_atoms = t.n_atoms
    n_dim = 3


    print('=== Calculating full (raw) msd data ===')
    msd_full_data = calc_msd_windowed(t,verbose=True)

    print('=== Evaluating Averages ===')
    msd_data = np.zeros([n_atoms,n_frames,7]) #to store summary for each chain/c.o.m.
    for ii in range(n_atoms):
        print('... working on c.o.m. {}'.format(ii))
        msd_data[ii,:] = calc_sliding_window_avg(msd_full_data[:,:,ii])

    print('=== Summarizing ===')
    msd_summary = np.mean(msd_data,0) #average over c.o.m. points
    msd_err = np.std(msd_data[:,:,1],0)/np.sqrt(n_atoms)
    msd_summary = np.insert(msd_summary,2,msd_err,axis=1)

    beta,beta_err = calc_D(msd_summary[:,0],msd_summary[:,1], err=msd_summary[:,2],n_dim=n_dim)

    print('=== Saving Results ===')
    np.savetxt('D_summary.txt', np.array([beta[1],beta_err[1], beta[0],beta_err[0]]), header='D\tD_err\ta\ta_err;\t\t msd = a + 2*ndim*Dt; ndim={}'.format(n_dim))

    np.savetxt('msd_summary.txt',msd_summary,header="msd\terr\twindow-std\twindow-err\twindow-t0\twindow-g\twindow-Neff;  msd = {}(+/-{}) + t*{}*{}(+/-{})".format(beta[0],beta_err[0],2*n_dim,beta[1],beta_err[1]))
    np.save('msd_full_data',msd_full_data)  #raw msd data at all time points
    np.save('msd_avg_data',msd_data)        #after sliding window averaging
    #note, the "window-X" quantities are chain-averaged values of all the sliding-window quantities. 
    #Get a sense for, on average, how each chain's own sliding window msd's are correlated

