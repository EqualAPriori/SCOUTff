import numpy as np

## RMSD calculation #DAN
def calc_rmsd(xyz,frame,atom_indices=None):
    if atom_indices is None:
        msds = np.sum((xyz - xyz[frame,:,:])**2.0,2) #subtract off frame 0, square, add along axis 2 (xyz coords)
    else:
        msds = np.sum((xyz - xyz[frame,atom_indices,:])**2.0,2) #subtract off frame 0, square, add along axis 2 (xyz coords)

    msd = np.mean(msds,1) #average across beads
    rmsd = np.sqrt(msd) # (nframe x 1) array, rmsd(t)
    return msd[frame:],rmsd[frame:]

## DAN
## need fitting function to extract diffusion from rmsd
# be careful to fit only to the linear part!
# dt = time step between successive frames
# delta_t = t-t_0, i.e. rmsd after a given period of time separation
# c x^2/(tau_ballistic + b*x)
# best is probably to plot and see what the curves look like
# ^^^ DO THAT, in case we get sub-diffusive regimes ^^^
def get_diffusion(rmsd,dt)

    return D


## decorrelation in diffusion times
# possibly use pymbar to remove the ends?
skip_last_num_frames = 10
Ds = np.zeros(traj.n_frames)
for iframe in xyzs.shape[0]-skip_last_num_frames:
    _,rmsd = calc_rmsd(xyz,iframe)
    Ds[iframe] = get_diffusion(rmsd) #gets an estimate of D using iframe as reference frame, i.e. D(t0)
# use pymbar to get decorrelation time of diffusion constant estimate



## block averaging
# in the below, need to be careful of indexing
# assume read in traj into n_frames
ndata = np.zeros(traj.n_frames-skip_last_num_frames)
msd_data = np.zeros([traj.n_frames,traj.n_frames-skip_last_num_frames]) #think of this as msd(t-t0)
n_frames = xyzs.shape[0]
for iframe in range(n_frames):
    msd,_ = calc_rmsd(xyz,iframe)
    msds[iframe,:msd.size] += msd
    ndata[:msd.size] += 1

for delta_t in range(n_frames):
    msd_tmp = msd_data[:ndata[delta_t],delta_t]
    #calculate statistics, i.e. mean, std, error in the mean; also get decorrelation time for each delta t!
    msd_mean[delta_t] =
    msd_std[delta_t] =
    msd_error[delta_t] = 

rmsd = np.sqrt(msd_mean)
#then fit diffusion constant


## NPT frame unwrapping #CHARLES









