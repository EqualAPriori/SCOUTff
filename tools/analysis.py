#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:03:04 2020

@author: My Nguyen
"""
import sys
sys.path.append('/home/mnguyen/bin/scripts/')
import scipy
import stats
import numpy as np
import matplotlib, sys, os
import matplotlib.pyplot as plt
import mdtraj as md
import log2txt
showPlots = True
try:
    os.environ["DISPLAY"] #Detects if display is available
except KeyError:
    showPlots = False
    matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window

N_av = 6.022e23 #1/mol
kB = 0.008314265 #kJ/mol/K

def GetThermo(ThermoLog, fi = 'lammps', obs = None, cols = None, autowarmup = True, warmup = 100, plot=False, plotDir = 'Thermo_plots'):
    """ fi: log file format, 'lammps' or 'openmm' """
    if not obs == None and not cols == None:
        Exception('Read data either by observable name or column index but not both!')

    if plot:
        try:
            os.mkdir(plotDir)
        except:
            pass
        print('...Thermo plots will be saved in {}...\n'.format(plotDir))
        
    #conver log file:
    if fi == 'openmm':
            ThermoLog = log2txt.log2txt_openmm([ThermoLog])[0]
    elif fi == 'lammps':
            section = 'PRODUCTION RUNS'
            ThermoLog = log2txt.log2txt_lammps([ThermoLog],section,'production')[0]

    print('new log file: {}'.format(ThermoLog))
    txt = ""
    obsID = []
    Stats = []
    #do stats
    file = open(ThermoLog,'r')
    if not obs == None:
        lines = file.readlines()
        while not isinstance(cols,list):
            for line in lines:
                if line.startswith('#'):
                                   obsNames = line.split()[1:]
                                   print('obsNames {}'.format(obsNames))
                                   cols = [obsNames.index(val) for val in obsNames if val in obs]
    print('cols {}'.format(cols))
    for i, col in enumerate(cols):
        if autowarmup:
            warmup,Data,nwarmup = stats.autoWarmupMSER(file, col)
            print ("Auto warmup detection with MSER-5 => ",nwarmup)
        else:
            warmup,Data = stats.extractData(file, col, warmup)
        (nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data, False ,False,'_{0}_mol{1}'.format(file.name,col))
        try:
            obsName = obsNames[col]
        except:
            obsName = 'col{}'.format(col)
        lines = "" 
        lines += '\n==== {} ===='.format(obsName)
        lines += "\n  - Mean                    = {} +/- {}".format(mean,semcc)
        lines += "\n  - Equilibrated samples    = {}".format(nsamples)
        lines += "\n  - Correlation time        = {}".format(kappa)
        lines += "\n  - Effective # samples     = {}".format(nsamples/kappa)
        lines += "\n  - Reduced-bias variance   = {}".format(unbiasedvar)
        # note that there is no unbiased estimator for the population standard deviation. We can use sqrt(var) as a indicative estimator.
        lines += "\n  - S.D. (unbiased, biased) = {} {}".format(np.sqrt(unbiasedvar),np.std(Data,ddof=0)) # ddof is correction to 1/N...using ddof=1 returns regular reduced-bias estimator
        lines += "\n  - Min, Max                = {} {}\n".format(min,max)

        ''' Plot '''
        if plot:
            plt.axvspan(0, nwarmup, alpha=0.5, color='#6495ED')
            plt.plot(np.hstack((warmup,Data)))
            plt.xlim(0)
            plt.xlabel('timestep')
            plt.ylabel(obsName)
            plt.savefig("{}/{}.png".format(plotDir,obsName),bbox_inches='tight')
            plt.close()
            
        print(lines)
        txt += lines
    
        Avg = mean
        Std = np.sqrt(unbiasedvar)
        Err = semcc
        CorrTime = kappa 
        NUncorrSamples = nsamples/kappa
        Stats.append([Avg,Std,CorrTime,Err,NUncorrSamples])
        obsID.append(obsName)

    return obsID, Stats
        
def GetRgRee(traj, DOP, NP, NAtomsPerChain = None, plotDir = 'RgRee_plots',
             RgDatName = 'RgTimeSeries', ReeDatName = 'ReeTimeSeries',RgStatOutName = 'RgReeStats', Ext='.dat',  
             res0Id = 0, autowarmup = True, warmup = 100, plot = False):    
    
    """NAtomsPerChain: used if running CG system, if provided will assume there is one residue per chain
       multiply coordinates by 10 if input traj was generated by lammps and unit is nonDim"""
    ElementDictionary ={
                        "carbon": 12.01,
                        "hydrogen": 1.008,
                        "oxygen": 16.00,
                        "nitrogen": 14.001,
                        "virtual site": 1.0,
                        "sodium": 23.0,
                        "chloride": 35.5}
    if plot:
        try:
            os.mkdir(plotDir)
        except:
            pass
        print('...Rg and Ree plots will be saved in {}...\n'.format(plotDir))
    RgTimeseries = [range(traj.n_frames)]
    Rgheader = "Frame   "
    
    RgSqStats = []
    RgSqTimeseries = [range(traj.n_frames)]
    RgSqheader = "Frame   "
    RgSqList = []

    txtRg = ""
    
    ReeTimeseries = [range(traj.n_frames)]
    Reeheader = "Frame   "
   
    ReeSqStats = []
    ReeSqTimeseries = [range(traj.n_frames)]
    ReeSqheader = "Frame   "
    ReeSqList = []

    #get indices of residues in all chains    
    MoleculeResidueList = []
    if not NAtomsPerChain:
        #number residues per chain = DOP (for AA systems)
        for j in range(NP):
            resId = range(res0Id + j*DOP, res0Id + (j+1)*DOP)
            MoleculeResidueList.append(resId)
    else:
        #1 residue per chain (for CG system)
        x = range(res0Id, res0Id + NP)
        MoleculeResidueList = [[a] for a in x]
        
    for j,resId in enumerate(MoleculeResidueList):
        resIdLow = np.min(resId)
        resIdUp = np.max(resId)
        atom_indices = traj.topology.select('resid {} to {}'.format(resIdLow,resIdUp)) 
        mass_list = []
        for index in atom_indices:
            element = str(traj.topology.atom(index).element)
            try:
                mass = ElementDictionary[element]
            except:
                mass = 1.
            mass_list.append(mass)
        mass_list = np.array(mass_list)
        if j == 0:
            print('Indices of atoms in chain {} \n{}'.format(j+1,atom_indices))       
            print('Mass of atoms in a chain {}'.format(mass_list))
        print('Evaluate Rg and Ree of chain {}/{}'.format(j+1,len(MoleculeResidueList)),end="\r")        

        '''=== Compute Rg ==='''
        Rg = md.compute_rg(traj.atom_slice(atom_indices),masses=mass_list) 
        RgTimeseries.append(Rg.tolist())
        Rgheader += 'Rg{}   '.format(j+1)
        np.savetxt(RgDatName+Ext, np.transpose(RgTimeseries), fmt = '%5.5f', header=Rgheader ) 

        RgSq = Rg**2.
        RgSqTimeseries.append(RgSq.tolist())  
        Rgheader += 'Rg{}^2   '.format(j+1)
        np.savetxt('RgSqTimeSeries'+Ext, np.transpose(RgSqTimeseries), fmt = '%5.5f', header=RgSqheader )
 
        #do stats on Rg^2
        file = open('RgSqTimeSeries'+Ext,'r')
        if autowarmup:
            warmup,Data,nwarmup = stats.autoWarmupMSER(file, j+1)
            #print ("Auto warmup detection with MSER-5 => ",nwarmup)
        else:
            nwarmup = warmup
            warmup,Data = stats.extractData(file, j+1, warmup)
        (nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data, False ,False,'_{0}_mol{1}'.format(file.name,j+1))
        Data = Data[::int(np.max([1.,kappa]))] # get decorrelated samples
        RgSqList.extend(Data)

        lines = ""
        lines += '\n==== Rg^2 for molecule {} ===='.format(j+1)
        lines += "\n  - Mean                    = {} +/- {}".format(mean,semcc)
        lines += "\n  - Equilibrated samples    = {}".format(nsamples)
        lines += "\n  - Correlation time        = {}".format(kappa)
        lines += "\n  - Effective # samples     = {}".format(nsamples/kappa)
        lines += "\n  - Reduced-bias variance   = {}".format(unbiasedvar)
        # note that there is no unbiased estimator for the population standard deviation. We can use sqrt(var) as a indicative estimator.
        lines += "\n  - S.D. (unbiased, biased) = {} {}".format(np.sqrt(unbiasedvar),np.std(Data,ddof=0)) # ddof is correction to 1/N...using ddof=1 returns regular reduced-bias estimator
        lines += "\n  - Min, Max                = {} {}\n".format(min,max)
        txtRg += lines

        Avg = mean
        Std = np.sqrt(unbiasedvar)
        Err = semcc
        CorrTime = kappa
        NUncorrSamples = nsamples/kappa
        RgSqStats.append([Avg,Std,CorrTime,Err,NUncorrSamples])

        ''' Plot Rg '''
        if plot:
            plt.axvspan(0, nwarmup, alpha=0.5, color='#6495ED')
            plt.plot(Rg, "k-")
            plt.xlim(0)
            plt.xlabel('timestep')
            plt.ylabel('Radius-of-gryation')
            plt.savefig("{}/Rg{}.png".format(plotDir,j+1),bbox_inches='tight')
            plt.close()

        '''=== Compute Ree ==='''
        atom_pairs = [np.min(atom_indices), np.max(atom_indices)]
        Ree = md.compute_distances(traj,atom_pairs= [atom_pairs], periodic=False, opt=True)
        Ree = Ree.tolist()
        Ree = [a[0] for a in Ree]
        ReeTimeseries.append(Ree)
        Reeheader += 'Ree{}   '.format(j+1)
        np.savetxt(ReeDatName+Ext, np.transpose(ReeTimeseries), fmt = '%5.5f', header=Reeheader ) 

        ReeSq = np.array(Ree)**2.
        ReeSqTimeseries.append(ReeSq.tolist())
        Reeheader += 'Ree{}^2   '.format(j+1)
        np.savetxt('ReeSqTimeSeries'+Ext, np.transpose(ReeSqTimeseries), fmt = '%5.5f', header=ReeSqheader )
               
        #do stats on Ree^2
        file = open('ReeSqTimeSeries'+Ext,'r')
        if autowarmup:
            warmup,Data,nwarmup = stats.autoWarmupMSER(file, j+1)
            #print ("Auto warmup detection with MSER-5 => ",nwarmup)
        else:
            nwarmup = warmup
            warmup,Data = stats.extractData(file, j+1, warmup)
        (nsamples,(min,max),mean,semcc,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data, False ,False,'_{0}_mol{1}'.format(file.name,j+1))
        Data = Data[::int(np.max([1.,kappa]))]
        ReeSqList.extend(Data)

        lines = ""
        lines += '\n==== Ree^2 for molecule {} ===='.format(j+1)
        lines += "\n  - Mean                    = {} +/- {}".format(mean,semcc)
        lines += "\n  - Equilibrated samples    = {}".format(nsamples)
        lines += "\n  - Correlation time        = {}".format(kappa)
        lines += "\n  - Effective # samples     = {}".format(nsamples/kappa)
        lines += "\n  - Reduced-bias variance   = {}".format(unbiasedvar)
        # note that there is no unbiased estimator for the population standard deviation. We can use sqrt(var) as a indicative estimator.
        lines += "\n  - S.D. (unbiased, biased) = {} {}".format(np.sqrt(unbiasedvar),np.std(Data,ddof=0)) # ddof is correction to 1/N...using ddof=1 returns regular reduced-bias estimator
        lines += "\n  - Min, Max                = {} {}\n".format(min,max)
        txtRg += lines

        Avg = mean
        Std = np.sqrt(unbiasedvar)
        Err = semcc
        CorrTime = kappa
        NUncorrSamples = nsamples/kappa
        ReeSqStats.append([Avg,Std,CorrTime,Err,NUncorrSamples])

        ''' Plot Ree '''
        if plot:
            plt.axvspan(0, nwarmup, alpha=0.5, color='#6495ED')
            plt.plot(Ree, "k-")
            plt.xlim(0)
            plt.xlabel('timestep')
            plt.ylabel('End-to-end distance')
            plt.savefig("{}/Ree{}.png".format(plotDir,j+1),bbox_inches='tight')
            plt.close()

    # get RMS Rg and Ree 
    RgSqList = np.array(RgSqList)
    RgRMS = np.sqrt(np.mean(RgSqList))
    RgSqErr = scipy.stats.sem(RgSqList)
    RgRMSErr = 1./2.*RgSqErr/RgRMS  # propagate SEM of Rg^2 to Rg
    RgSqStd = np.std(RgSqList,ddof=1)
    RgRMSStd = 1./2.*RgSqStd/RgRMS  # propagate Std of Rg^2 to Rg 
    RgSqStats = np.array(RgSqStats)
    RgRMSCorrTime = np.mean(RgSqStats[:,2])
    RgRMSCorrTimeErr = np.sqrt(np.var(RgSqStats[:,2])/len(RgSqStats[:,2]))
    RgRMSNUncorrSamples = np.mean(RgSqStats[:,4])
     
    ReeSqList = np.array(ReeSqList)
    ReeRMS = np.sqrt(np.mean(ReeSqList))
    ReeSqErr = scipy.stats.sem(ReeSqList)
    ReeRMSErr = 1./2.*ReeSqErr/ReeRMS 
    ReeSqStd = np.std(ReeSqList,ddof=1)
    ReeRMSStd = 1./2.*ReeSqStd/ReeRMS  
    ReeSqStats = np.array(ReeSqStats)
    ReeRMSCorrTime = np.mean(ReeSqStats[:,2])
    ReeRMSCorrTimeErr = np.sqrt(np.var(ReeSqStats[:,2])/len(ReeSqStats[:,2]))
    ReeRMSNUncorrSamples = np.mean(ReeSqStats[:,4])
    
    lines = ""
    lines += '\n\n====================='
    lines += '\n\nRMS of Rg is: {0:2.4f} +/- {1:2.5f}'.format(RgRMS, RgRMSErr)
    lines += '\nRMS Rg correlation time: {0:5.4f} +/- {1:5.6f}'.format(RgRMSCorrTime, RgRMSCorrTimeErr)
    lines += '\n\nRMS of Ree is: {0:2.4f} +/- {1:2.5f}'.format(ReeRMS, ReeRMSErr)
    lines += '\nRMS Ree correlation time: {0:5.4f} +/- {1:5.6f}'.format(ReeRMSCorrTime, ReeRMSCorrTimeErr)
    
    print(lines)
    txtRg += lines
    f = open(RgStatOutName+Ext,'w')
    f.write(txtRg)
    return RgRMS,ReeRMS,RgRMSErr,ReeRMSErr,RgRMSCorrTime,RgRMSCorrTimeErr,RgRMSNUncorrSamples,ReeRMSCorrTime,ReeRMSCorrTimeErr,ReeRMSNUncorrSamples,RgRMSStd,ReeRMSStd

def RMSBond(traj, backboneAtoms = ['C1','C2'], resid = [], autowarmup=True, warmup=100):
    """backboneAtoms: atom names of heavy backbone atoms in the topology"""
    
    top = traj.topology
    
    # slice trajectory if provide resid, otherwise calculate bond length on the full trajectory
    if len(resid)>0:
        aId = []
        for i in resid:
            aId.extend(top.select("resid {}".format(i)))
        traj = traj.atom_slice(aId)
        top = traj.topology
        
    # get pairs of atoms to calculate bond length
    IsBackbone = [set([bond[0].name,bond[1].name]).issubset(set(backboneAtoms)) for bond in top.bonds]
    IsBackbone = np.array(IsBackbone,bool)
    pairs = np.array([[bond[0].index,bond[1].index] for bond in top.bonds])
    pairs = pairs[IsBackbone]
    
    l = md.compute_distances(traj,pairs)
    l = l.flatten()
        
    lAvg = np.mean(l)
    lStd = np.std(l)
    lErr = scipy.stats.sem(l)
    
    # get RMS bond length     
    lRMS = np.sqrt(np.mean(l**2.))
    lSqErr = scipy.stats.sem(l**2)
    lRMSErr = 1./2.*lSqErr/lRMS
    
    print('\n\nRMS bond length: {0:2.5f} +/- {1:2.5f}'.format(lRMS,lRMSErr))
    print('from heavy atoms {}\n'.format(backboneAtoms))
    return lAvg,lStd,lErr, lRMS, lRMSErr     
        
def GetStats(trajFile, top, NP, ThermoLog, DOP = 10, NAtomsPerChain = None,  
             backboneAtoms= ['C1','C2'], monMass={}, nMons={}, density=None,
             StatsFName = 'AllStats.dat', RgDatName = 'RgTimeSeries', ReeDatName = 'ReeTimeSeries',RgStatOutName = 'RgReeStats', Ext='.dat',  
             fi = 'openmm', obs = None, cols = None, density_col = None,
             res0Id = 0, stride = 1, autowarmup = True, warmup = 100, plot = False, unit = 'real'):
    """"
    fi = 'openmm' or 'lammps'
    """
    txt = '# Obs.    Avg.\tS.D.\tStdErr.\tCorr.\tStdErr.\tUncorr.Samples\n'
    traj = md.load(trajFile, top=top, stride = stride)
    traj.make_molecules_whole(inplace=True, sorted_bonds=None)
    if fi == 'lammps' and unit == 'nonDim':
        traj.xyz *= 10.
        traj.unitcell_lengths *= 10
    
    if NP > 0:
        resid = range(res0Id,res0Id +  NP * DOP - 1) # residues indices of polymer chains
        RgRMS,ReeRMS,RgRMSErr,ReeRMSErr,RgRMSCorrTime,RgRMSCorrTimeErr,RgRMSNUncorrSamples,ReeRMSCorrTime,ReeRMSCorrTimeErr,ReeRMSNUncorrSamples,RgRMSStd,ReeRMSStd = GetRgRee(traj, DOP, NP, NAtomsPerChain = NAtomsPerChain,
             RgDatName = RgDatName, ReeDatName = ReeDatName, RgStatOutName = RgStatOutName, Ext=Ext,
             res0Id = res0Id, autowarmup = autowarmup, warmup = warmup, plot = plot)

        txt += 'RMSRg  %8.5f  %8.5f  %8.5f %8.5f  %8.5f  %i'%(RgRMS,RgRMSStd,RgRMSErr,RgRMSCorrTime,RgRMSCorrTimeErr,RgRMSNUncorrSamples)
        txt += '\nRMSRee  %8.5f  %8.5f  %8.5f %8.5f  %8.5f  %i'%(ReeRMS,ReeRMSStd,ReeRMSErr,ReeRMSCorrTime,ReeRMSCorrTimeErr,ReeRMSNUncorrSamples)
        
    if density_col != None:
        if cols == None:
            cols = []
        cols.append(density_col)

    if not (obs==None and cols==None):
        print('reading thermo file {}'.format(ThermoLog))
        obsID, Stats = GetThermo(ThermoLog, fi = fi, obs = obs, cols = cols, autowarmup = autowarmup, warmup = warmup, plot=plot)    
        for i, obs in enumerate(obsID):
            Avg,Std,CorrTime,Err,NUncorrSamples = Stats[i]
            if density_col != None and i == len(cols)-1:
                obsName = 'Density'
                density = Avg
            else: 
                obsName = obs
            try:
                txt +=  '\n%s  %8.5f  %8.5f  %8.5f  %8.5f  %s  %i' %(obsName, Avg, Std, Err, CorrTime, 'N/A',NUncorrSamples)
            except:
                txt +=  '\n%s  %8.5f  %8.5f  %8.5f  %8.5f  %s  %s' %(obsName, Avg, Std, Err, CorrTime, 'N/A',NUncorrSamples)
    txt2 = ""
    if backboneAtoms and NP > 0:
        lAvg,lStd,lErr,lRMS, lRMSErr = RMSBond(traj, backboneAtoms=backboneAtoms, resid=resid)
        txt +=  '\n%s  %8.8f  %8.8f  %8.5f' %('Bond', lAvg, lStd, lErr)
        #calculate characteristic ratio
        n = len(backboneAtoms)*DOP - 1 # number of backbone bonds
        Cn = ReeRMS**2 / (n * lRMS**2)
        txt2 += '\n\nRMS bond length between heavy atoms: %8.8f +/- %8.8f' %(lRMS,lRMSErr)
        txt2 += '\nCharacteristic ratio:  %8.8f' %(Cn)
        
    # statistical segment length
    Vref1 = 0.1 #nm^3
    Vref2 = 0.117 #nm^3
    if density and NP > 0 and len(nMons.items())>0:
        #check if DOP = total number of monomers:
        if np.abs(np.sum(np.array(list(nMons.values())))-DOP) > 1e-2: 
            raise Exception('Total number of monomers in monInChain must be equal to DOP')
        chainMass = 0 # g/mol
        for mon, n in nMons.items():
            chainMass += float(n) * float(monMass[mon])
        chainMass += 2. # 2 extra H's for the ends
        chainMass /= N_av # g/chain
        Vchain = chainMass/density * 1.0e21 # nm^3/chain
        Nseg1 = Vchain/Vref1
        Nseg2 = Vchain/Vref2
        b1 = RgRMS * np.sqrt(6./Nseg1)
        b2 = RgRMS * np.sqrt(6./Nseg2)
        txt2 += '\nb (Vref=%3.4fnm^3): %8.8f' %(Vref1,b1)
        txt2 += '\nb (Vref=%3.4fnm^3): %8.8f' %(Vref2,b2)
    print(txt2)
    f = open(StatsFName, 'w')
    print('\n...Writing results to {}...'.format(StatsFName))
    f.write(txt)
    f.write(txt2)
    f.close()

def GetCompressibility(trajFile, top, temp, stride = 1, unit = 'bar', lengthScale = 0.31, fName = 'Compressibility', Ext='.dat', trajFmt = 'omm'):
    """unit = ['bar','Pa','nonDim']
       lengthScale in nm
       trajFmt: lmp or omm, multiply volume by 10**3 if trajFmt is lmp and unit is nonDim"""

    print('Need temperature in Kelvin')
    kT = kB * temp #kJ/mol 
    traj = md.load(trajFile, top=top, stride = stride)
    vols = traj.unitcell_volumes
    if trajFmt == 'lmp' and unit =='nonDim':
        vols *= 10.**3
    meanVol = np.mean(vols)
    vol2s = vols**2
    meanVol2 = np.mean(vol2s)
    if unit == 'Pa':
        compressibility = (meanVol2 - meanVol**2)/(meanVol*kT) * N_av * 1e-30
        s = '1/'+unit
        print('Compressibility is {:.4e} {:s}'.format(compressibility, s))
    elif unit == 'bar':
        compressibility = (meanVol2 - meanVol**2)/(meanVol*kT) * N_av * 1e-25
        s = '1/'+unit
        print('Compressibility is {:.4e} {:s}'.format(compressibility, s)) 
    elif unit == 'nonDim':
        kT = 1.
        print('Assume kT = 1.')
        compressibility = (meanVol2 - meanVol**2)/(meanVol*kT) * N_av
        compressibilityBar = compressibility * lengthScale**3 / (kB * temp) * 1e-25
        s = 'sigma^3/kT'
        print('Compressibility is {:.4e} {:s}'.format(compressibility, s))
        print('Convert to real unit at {:.2f} K for length scale of {:.2f} nm: {:.4e} 1/bar'.format(temp, lengthScale, compressibilityBar))
    f = open(fName+Ext,'w')
    f.write('{:.5e} {:s}'.format(compressibility, s))
    return compressibility

if __name__ ==  '__main__':
    
    import argparse as ap
    parser = ap.ArgumentParser()
    parser.add_argument('-traj', required=True, type=str, help = "trajectory file")
    parser.add_argument('-top', required=True, type=str, help = "topology file")
    parser.add_argument('-s', type=int, default=1, help = "trajectory stride")
    parser.add_argument('-ther', type=str, default=None, help = 'thermo log file')
    parser.add_argument('-np', type = int, required = True, default = 1, help = 'number of polymer chains')
    parser.add_argument('-dop', type = int, required = True, default = 1, help = 'chain length')
    parser.add_argument('-c', type = int, nargs ='+', default = None, help = 'column indices (index of 1st column is 0) in thermo log file to get stats')
    parser.add_argument('-densc', type = int, default = None, help = 'column indices of melt density in thermo log file, must be in g/mL, override -dens if both are  provided')
    parser.add_argument('-dens', type = float, default = None, help = 'melt density in g/mL')
    parser.add_argument('-bb', type = str, default = None, nargs ='+', help='names of backbone atoms in the topology file')                            
    parser.add_argument('-mon', action='append', nargs=2, help = "pairs: monomer_name number_per_chain")

    parser.add_argument('-res0', type = int, default = 0, help='index of the first polymer residue')
    parser.add_argument('-a', action='store_true', help='use autowarmup to get stats')
    parser.add_argument('-g', action='store_true', help='plot time series of observables')
    parser.add_argument('-w', type = int, default=100, help='number of warmup data points')    
    args = parser.parse_args()
    
    """ example command:
    python analysis.py -traj run0_output.dcd -top run0_post_production.pdb -ther run0_thermo_production.out -np 50 -dop 60 -c 4 -densc 6 -bb C1 C2 -mon s12pB 60 -a -g"""

    TrajFile = args.traj  
    top = args.top  
    ThermoLog = args.ther    
    stride = args.s
    
    # Rg Ree    
    NP = args.np
    DOP = args.dop 
    res0Id = args.res0
    
    # thermo 
    cols = args.c
    density_col = args.densc
    density = args.dens
    
    # RMS bond between adjacent heavy atoms in the backbone
    backboneAtoms = args.bb
    
    # segment length from density (when density_col or density is provided)
    # list monomers in the chain and their respective numbers (monomers must have mass defined in monomerMass dictionary):
    nMons = {}
    if args.mon:
        for val in args.mon:
            nMons.update({val[0]:int(val[1])})
    monMass = {'pS': 108., 'pmS': 118., 's12pB': 56., 's14pB': 56., '12pB': 54., '14pBcis': 54., '14pBtrans':54.} # g/mol
    
    autowarmup = args.a
    plot = args.g
    warmup = args.w
    
    GetStats(TrajFile, top, NP, ThermoLog, DOP = DOP, NAtomsPerChain = None, stride = stride, 
            backboneAtoms = backboneAtoms, StatsFName = 'AllStats.dat', 
            monMass=monMass,nMons=nMons, density=density, res0Id = res0Id,
            RgDatName = 'RgTimeSeries', ReeDatName = 'ReeTimeSeries',RgStatOutName = 'RgReeStats', Ext='.dat',
            cols = cols, density_col = density_col, autowarmup = autowarmup, warmup = warmup, plot = plot)

