#! Helper functions for alchemical calculation 
"""
Updated  2021.01.13
@author: Kevin Shen
"""

""" Summary of components:
1) function and/or object for modifying forces in the system
2) function for recalculating energies
3) function for assembling data and calculating free energy difference
4) also will need to write a separatte script for doing the file i/o and calling all these functions
"""

import mdtraj
from pymbar import MBAR, timeseries
import simtk.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as unit


def add_cnb(system,indices,suffix='',lambda0=1.0):
    """
    takes in a system, turns off alchemical
    only adds a single lambda controlling parameter!

    suffix is for adding a modifier to change the name of the lambda parameter
    """

    # Retrieve the NonbondedForce
    forces = { force.__class__.__name__ : force for force in system.getForces() }
    nbforce = forces['NonbondedForce']

    # Add a CustomNonbondedForce to handle only alchemically-modified interactions
    alchemical_particles = set(list(indices))
    chemical_particles = set(range(system.getNumParticles())) - alchemical_particles
    energy_function = 'lambda{}*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;'.format(suffix)
    energy_function += 'reff_sterics = sigma*(0.5*(1.0-lambda{}) + (r/sigma)^6)^(1/6);'.format(suffix)
    energy_function += 'sigma = 0.5*(sigma1+sigma2); epsilon = sqrt(epsilon1*epsilon2);'
    custom_force = openmm.CustomNonbondedForce(energy_function)
    custom_force.addGlobalParameter('lambda{}'.format(suffix), lambda0)
    custom_force.addPerParticleParameter('sigma')
    custom_force.addPerParticleParameter('epsilon')
    for index in range(system.getNumParticles()):
        [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
        custom_force.addParticle([sigma, epsilon])
        if index in alchemical_particles:
            nbforce.setParticleParameters(index, charge*0, sigma, epsilon*0)
     
    custom_force.addInteractionGroup(alchemical_particles, chemical_particles)
    system.addForce(custom_force)
#end def add_cnb()


def update_state(context, states):
    """
    states should be a list of 2-tuples of (lambda_name, lambda_value)
    """
    for state in states:
        lambda_name, lambda_value = state
        context.setParameter(lambda_name, lambda_value)


def recalc_energy(context, traj):
    """
    initialize array for storing energies
    traj is the reference trajectory
    """
    PE = np.zeros(traj.n_frames) 
    for iframe,frame in enumerate(traj):
        if np.mod(iframe,100) == 0:
            logger.info("...Frame {}".format(it))
        box = frame.unitcell_vectors[0]
        context.setPeriodicBoxVectors(box[0], box[1], box[2])
        context.setPositions(frame.xyz[0])
        state = context.getState(getEnergy=True)
        PE[iframe] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        
    with open(outfile,'w') as f:
        f.write("#frame\tPE(kJ/mol), ewald error tolerance: {}\n".format(args.ewald_error_tolerance))
        for ie,energy in enumerate(PE):
            f.write("{}\t{}\n".format(ie,energy))


def calc_df(u_kln):
    """
    u_kln should be (nstates) x (nstates) x (nframes)
    note that u_kln should be normalized by kT already
    where each element is 
        a config from frame `n` of a trajectory conducted with state `k`
        with energy recalculated using parameters of state `l`
    """
    dims = u_kln.shape
    if dims[0] != dims[1]:
        raise ValueError("dimensions {} of u_kln should be square in the first two indices".format(dims))
    nstates = dims[0]

    N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
    for k in range(nstates):
        [nequil, g, Neff_max] = timeseries.detectEquilibration(u_kln[k,k,:])
        indices = timeseries.subsampleCorrelatedData(u_kln[k,k,:], g=g)
        N_k[k] = len(indices)
        u_kln[k,:,0:N_k[k]] = u_kln[k,:,indices].T
    # Compute free energy differences and statistical uncertainties
    mbar = MBAR(u_kln, N_k)
    [DeltaF_ij, dDeltaF_ij, Theta_ij] = mbar.getFreeEnergyDifferences()

    # save data?

    return DeltaF_ij, dDeltaF_ij


def assemble_df(u_kln,mode=0):
    """
    Should suppor two modes:
    1) have the full u_kln matrix of eneries from state k to state l
    2) only have energy differences from adjacent states (i.e. all the other elements, even if missing, are ignored)

    Note that this function right now assumes data is already assembled into a u_kln matrix... so will also need helper script or function to do that
    """
    if mode == 0:
        print('=== mode 0, using full u_kln matrix ===')
        dF, ddF = df(u_kln)
    else:
        nstates = u_kln.shape[0]
        dFs = np.zeros(nstates-1)
        ddFs = np.zeros(nstates-1)
        for state in range(nstates-1):
            u_kln_tmp = u_kln[state:state+2, state:state+2, :]
            dFs_tmp, ddFs_tmp = df(u_kln)
            dFs[state] = Fs_tmp[0,1]
            ddFs[state] = ddFs_tmp[0,1]
            
        np.savetxt('DeltaF_{}.dat'.format(flag),dFs)
        np.savetxt('dDeltaF_{}.dat'.format(flag),ddFs)

        # Print out one line summary 
        with open('muex_{}.dat'.format(flag),'w') as f:
                f.write('DeltaF:\t {} +/- {} kT'.format( dFs.sum(), np.sqrt(np.sum(ddFs**2)) ) )
                print('DeltaF:\t {} +/- {} kT'.format( dFs.sum(), np.sqrt(np.sum(ddFs**2)) ) )



