import os
import sys
import time
import numpy as np
import shutil
import mdtraj as md

''' ProgressBar Class '''

ProgressBarOn = sys.stdout.isatty()

class ProgressBar(object):
    def __init__(self, Text, Steps = 1, BarLen = 20, UpdateFreq = 1.):
        """Initializes a generic progress bar."""
        self.Text = Text
        self.Steps = Steps
        self.BarLen = BarLen
        self.UpdateFreq = UpdateFreq
        self.__LastTime = 0.
        self.__LastTime = time.time()
        self.__LastLen = 0
        self.Update(0)

    def Update(self, Step):
        """Updates the progress bar."""
        if time.time() - self.__LastTime > self.UpdateFreq:
            if not ProgressBarOn:
                return
            self.__LastTime = time.time()
            if self.BarLen == 0:
                s = "%s [%d]" % (self.Text, Step)
            else:
                Frac = float(Step) / (self.Steps + 1.e-300)
                n = int(self.BarLen * Frac + 0.5)
                n = max(min(n, self.BarLen), 0)
                s = "%s [" % self.Text
                s += "="*n + (self.BarLen-n)*" "
                s += "] %.2f%%" % (100.*Frac)
            self.__LastLen = len(s)
            s += "\r"
            sys.stdout.write(s)
            sys.stdout.flush()

    def Clear(self):
        """Clears text on this line."""
        if not ProgressBarOn:
            return
        sys.stdout.write(" "*self.__LastLen + "\r")
        sys.stdout.flush()

class TrajObj():
    ''' Object that loads in a trajectory 
    
    
    
    
    
    '''
    
    
    def __init__(self,_TrajList,_Topology,_Warmup=0,_Stride=0.,
                    _CreateSaveDir=False,_SaveDirName='TrajMap'):
        ''' 
            Input:
                
                (1) _Topology - topology file for trajectories
                (2) _TrajList - a list of trajectories to load
                (3) _Warmup   - number of frames to skip
                (4) _Stride   - traj load frequency
                (5) _CreateSaveDir - Save outputs to a file
        
        '''
        
        
        self.TrajList       = _TrajList 
        self.Topology       = _Topology
        #self.nSpecies       = int(_nSpecies) currently not used
        self.Warmup         = int(_Warmup)
        self.Stride         = int(_Stride) 
        self.Verbose        = False
        self.RemoveSolvent  = False
        
        self.ResidueNames   = [] # list of unique residue names
        self.LogFileName    = 'TrajMap.Log'
        self.SaveDirName    = _SaveDirName
        self.CreateSaveDir  = bool(_CreateSaveDir)
        self.SavePath       = os.path.join(os.getcwd(),self.SaveDirName)    
        self.SetSaveDir(self.SaveDirName,self.CreateSaveDir)    
        self.ResidueMasses  = '' # dictionary of np.arrays of atom masses for each residue
        self.WrapTraj       = True
        self.nProcessors    = int(1)     
        
        # overwrite old log file
        open(self.LogFileName,'w')
        
        # load in trajectory
        self.mdTraj     = md.load(self.TrajList,top=self.Topology)
        if self.RemoveSolvent: # for debugging, etc.
            self.mdTraj.remove_solvent(inplace=True)
        self.mdTraj     = self.mdTraj[self.Warmup:]
        self.mdTraj     = self.mdTraj[::self.Stride]
        self.Write2Log('Loaded Trajectories:\n')
        self.Write2Log('{}\n'.format(self.TrajList))
        self.Write2Log('Traj. Warmup:    {}\n'.format(self.Warmup))
        self.Write2Log('Traj. Slice:     {}\n'.format(self.Stride))
        self.Write2Log('Number Frames:   {}\n'.format(self.mdTraj.n_frames))
        self.Write2Log('Number Chains:   {}\n'.format(self.mdTraj.n_chains))
        self.Write2Log('Number Residues: {}\n'.format(self.mdTraj.n_residues))
        self.Write2Log('Number Atoms:    {}\n'.format(self.mdTraj.n_atoms))
        
        self.getResidues()
        self.Write2Log('\n')
        self.Write2Log('Residue Masses:\n')
        for _indx,resname in enumerate(self.ResidueNames):
            self.Write2Log('Residue {}:   {}\n'.format(resname,self.ResidueMasses[resname]))
            self.Write2Log('Mass Sum: {}\n'.format(np.sum(self.ResidueMasses[resname])))
            self.Write2Log('\n')
        
        
    def SetSaveDir(self,_SaveDirName,_CreateSave=False):
        ''' Set Save Directory Name '''
        if _CreateSave:
            self.SaveDirName = str(_SaveDirName)
            try:
                os.mkdir(self.SaveDirName)
            except: # overwriting old directory
                shutil.rmtree(self.SaveDirName)
                os.mkdir(self.SaveDirName)
            self.LogFileName = os.path.join(os.getcwd(),self.SaveDirName,self.LogFileName)
            
    def getResidues(self):
        ''' 
            Builds list of the names of the residues, and 
            collects masses of atoms for each residue. 
        '''
        
        _temp = [] # residue list
        _mass_temp = {} # residue dictionary of mass arrays
        for _indx,residue in enumerate(self.mdTraj.topology.residues):
            if residue.name in _temp:
                pass
            else:
                _temp.append(residue.name)
                # collect masses for each atom in the residue
                _mtemp = []
                for _aindx,atom in enumerate(residue.atoms):
                    _mtemp.append(atom.element.mass)
                
                _mass_temp['{}'.format(residue.name)] = np.asarray(_mtemp)
        
        self.Write2Log('\n')
        self.Write2Log('Finding Unique Residues:\n')        
        self.Write2Log('{}\n'.format(_temp))
        self.ResidueNames = _temp
        self.ResidueMasses = _mass_temp
   
    def Write2Log(self,_text):
        """ Write text to log file """
        try:
            _f = open(self.LogFileName,'a')
        except:
            _f = open(self.LogFileName,'w')
        
        _f.write('{}'.format(_text))
        _f.close()
        
    def DoCOMMap(self):
        ''' Function to perform COM mapping of Loaded Trajectories 
                based off of residues.

        '''
        
        PBar = ProgressBar('Mapping Traj:', Steps = (int(self.mdTraj.n_residues)-1), BarLen = 20, UpdateFreq = 1.)
        
        self.Write2Log('\n')
        self.Write2Log('Performing Mapping:\n')
       
       # setup topology file for COM trajectory
        _top2 = md.Topology()
        for rindx, residue in enumerate(self.mdTraj.topology.residues):
            _top2.add_chain()
            _top2.add_residue(residue.name,_top2.chain(-1))
            _top2.add_atom("P",md.element.carbon,_top2.residue(-1))
        
        start = time.time()
        # setup np matrix to hold the mapped coordinates
        COM_xyz = np.zeros((int(self.mdTraj.n_frames),int(self.mdTraj.n_residues),3))
        # loop through each residue and do mapping for ALL frames in trajectory
        for rindx, residue in enumerate(self.mdTraj.topology.residues):
            PBar.Update(rindx)
            resTraj = self.mdTraj.atom_slice(self.mdTraj.topology.select('resid {}'.format(rindx)))
            tempmass = np.column_stack((self.ResidueMasses[residue.name],self.ResidueMasses[residue.name],self.ResidueMasses[residue.name]))
            rsq = np.multiply(resTraj.xyz,tempmass)
            rsq = np.sum(rsq,axis=1)/np.sum(self.ResidueMasses[residue.name])
            if self.WrapTraj: # WrapTrajCoordinates
                rsq = np.mod(rsq,self.mdTraj.unitcell_lengths)
            
            COM_xyz[:,rindx,:] = rsq
        
        final = time.time()
        totaltime = final - start
        self.Write2Log('RunTime:            {}\n'.format(totaltime))
        PBar.Clear()
        
        trajCOM = md.Trajectory(COM_xyz,_top2,unitcell_lengths = self.mdTraj.unitcell_lengths, unitcell_angles = self.mdTraj.unitcell_angles)
        trajCOM.save_lammpstrj(os.path.join(self.SaveDirName,"TrajCOM.lammpstrj"))
        trajCOM.save_dcd(os.path.join(self.SaveDirName,"TrajCOM.dcd"))
        _PDB = md.formats.PDBTrajectoryFile(os.path.join(self.SaveDirName,"TrajCOM.pdb"), mode='w', force_overwrite=True, standard_names=True)
        _PDB.write(trajCOM.xyz[0],trajCOM.topology)
        
    def DoCOMMapPLL_OLD(self):
        ''' Function to perform COM mapping of Loaded Trajectories 
                based off of residues in parallel.
            
            Not Used!
        '''
        import threading
        import logging
        import multiprocessing
        from joblib import Parallel, delayed
        
        PBar = ProgressBar('Mapping Traj:', Steps = (int(self.mdTraj.n_residues)-1), BarLen = 20, UpdateFreq = 1.)
        
        # setup topology file for COM trajectory
        _top2 = md.Topology()
        for rindx, residue in enumerate(self.mdTraj.topology.residues):
            _top2.add_chain()
            _top2.add_residue(residue.name,_top2.chain(-1))
            _top2.add_atom("P",md.element.carbon,_top2.residue(-1))
        
        # setup np matrix to hold the mapped coordinates
        COM_xyz = np.zeros((int(self.mdTraj.n_frames),int(self.mdTraj.n_residues),3))
        _cnt = 0
        def DoCOMMap(_range,mdTraj,ResidueMasses,COM_xyz,cnt):

            for rindx in _range:
                resTraj = mdTraj.atom_slice(mdTraj.topology.select('resid {}'.format(rindx)))
                residue = mdTraj.topology.residue(rindx)
                tempmass = np.column_stack((ResidueMasses[residue.name],ResidueMasses[residue.name],self.ResidueMasses[residue.name]))
                rsq = np.multiply(resTraj.xyz,tempmass)
                rsq = np.sum(rsq,axis=1)/np.sum(ResidueMasses[residue.name])
                COM_xyz[:,rindx,:] = rsq
                cnt += 1
                PBar.Update(cnt)
        
        
        _nThreads = 4
        _intperthread = int(self.mdTraj.n_residues/_nThreads)
        _remainder = np.mod(self.mdTraj.n_residues,_nThreads)
        _nPerThread = []
        for _i in range(_nThreads):
            if _i != _nThreads-1:
                _nPerThread.append(_intperthread)
            else:
                _nPerThread.append(_intperthread+_remainder)
        
        
        self.Write2Log('\n')
        self.Write2Log('Threading Report:\n')
        self.Write2Log('Total Number Residues: {}\n'.format(sum(_nPerThread)))
        self.Write2Log('nPerThread:            {}\n'.format(_nPerThread))
        
        start = time.time()
        threads = []        
        for _i in range(_nThreads):
            _npt = _nPerThread[_i]
            _range = np.linspace((_i*_npt),(_i+1)*_npt-1,_npt,dtype='int64')
            print('RANGE')
            print(_range)
            _t = threading.Thread(target=DoCOMMap, args=(_range,self.mdTraj,self.ResidueMasses,COM_xyz,_cnt))
            threads.append(_t)
            _t.start()
        
        for index, thread in enumerate(threads):
            logging.info("Main    : before joining thread %d.", index)
            thread.join()
            logging.info("Main    : thread %d done", index)
        
        PBar.Clear()
        final = time.time()
        totaltime = final - start
        self.Write2Log('RunTime:            {}\n'.format(totaltime))
        
        
        trajCOM = md.Trajectory(COM_xyz,_top2,unitcell_lengths = self.mdTraj.unitcell_lengths, unitcell_angles = self.mdTraj.unitcell_angles)
        trajCOM.save_lammpstrj(os.path.join(self.SaveDirName,"TrajCOM.lammpstrj"))
        trajCOM.save_dcd(os.path.join(self.SaveDirName,"TrajCOM.dcd"))
        _PDB = md.formats.PDBTrajectoryFile(os.path.join(self.SaveDirName,"TrajCOM.pdb"), mode='w', force_overwrite=True, standard_names=True)
        _PDB.write(trajCOM.xyz[0],trajCOM.topology)
                
    def DoCOMMapPLL(self,_nProcessors=1):
        ''' Function to perform COM mapping of Loaded Trajectories 
                based off of residues in parallel.

        '''
        import threading
        import logging
        import multiprocessing as mp
        import ctypes
        lock = mp.Lock()
        #manager = mp.Manager()
        #man_list = manager.list()
        
        self.nProcessors = int(_nProcessors)
        PBar = ProgressBar('Mapping Traj:', Steps = (int(self.mdTraj.n_residues)-1), BarLen = 20, UpdateFreq = 1.)
        
        # setup topology file for COM trajectory
        _top2 = md.Topology()
        for rindx, residue in enumerate(self.mdTraj.topology.residues):
            _top2.add_chain()
            _top2.add_residue(residue.name,_top2.chain(-1))
            _top2.add_atom("P",md.element.carbon,_top2.residue(-1))
        
        # setup np matrix to hold the mapped coordinates
        COM_xyz = np.zeros((int(self.mdTraj.n_frames),int(self.mdTraj.n_residues),3))
        _cnt = 0
        def DoCOMMap(_range,mdTraj,ResidueMasses,COM_xyz,_cnt,lock):
            for rindx in _range:
                try:
                    resTraj = mdTraj.atom_slice(mdTraj.topology.select('resid {}'.format(rindx)))
                except:
                    print('Error Slice Residue: {}...Skipping'.format(rindx))
                    pass
                residue = mdTraj.topology.residue(rindx)
                tempmass = np.column_stack((ResidueMasses[residue.name],ResidueMasses[residue.name],ResidueMasses[residue.name]))
                rsq = np.multiply(resTraj.xyz,tempmass)
                rsq = np.sum(rsq,axis=1)/np.sum(ResidueMasses[residue.name])
                if self.WrapTraj: # WrapTrajCoordinates
                    rsq = np.mod(rsq,self.mdTraj.unitcell_lengths)
                lock.acquire() # do not let other process write to same array
                COM_xyz[:,rindx,:] = rsq
                _cnt.value += 1
                PBar.Update(_cnt_share.value)
                lock.release() # release this processes hold
                
        
        # get the number of residues per process
        _nThreads = self.nProcessors
        _intperthread = int(self.mdTraj.n_residues/_nThreads)
        _remainder = np.mod(self.mdTraj.n_residues,_nThreads)
        _nPerThread = []
        for _i in range(_nThreads):
            if _i != _nThreads-1:
                _nPerThread.append(_intperthread)
            else:
                _nPerThread.append(_intperthread+_remainder)
        
        # output some stuff about the multiprocess job
        self.Write2Log('\n')
        self.Write2Log('Threading Report:\n')
        num_cores = mp.cpu_count()
        self.Write2Log('Number Cores:          {}\n'.format(num_cores))
        self.Write2Log('Total Number Residues: {}\n'.format(sum(_nPerThread)))
        self.Write2Log('nPerThread:            {}\n'.format(_nPerThread))
        
        start = time.time()
        processes = []  
        temp_range = []
        
        # Create Shared Memory objects
        _cnt_share = mp.Value('i',0) # shared counter
        
        # Nice little function for helping setup shared arrays that work with numpy, maybe not the cleanest
        def shared_array(shape):
            """
            Form a shared memory numpy array.
            
            http://stackoverflow.com/questions/5549190/is-shared-readonly-data-copied-to-different-processes-for-python-multiprocessing 
            """
    
            shared_array_base = mp.Array(ctypes.c_double, shape[0]*shape[1]*shape[2])
            shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
            shared_array = shared_array.reshape(*shape)
            return shared_array


        # Form a shared array and a lock, to protect access to shared memory.
        array = shared_array(COM_xyz.shape) # shared trajectory positions
        
        # Split up the residues onto the different processes 
        _npt_current = 0
        for _i in range(_nThreads):
            _npt = _nPerThread[_i] 
            _range = np.linspace(_npt_current,_npt_current+_npt-1,_npt,dtype='int64')
            _npt_current += _npt
            if self.Verbose: # for troubleshooting
                self.Write2Log("\n")
                self.Write2Log('Range of Residues for Process {}\n'.format(_i))
                self.Write2Log("{}\n".format(_range))
            temp_range.append(_range)
            _p = mp.Process(target=DoCOMMap, args=(_range,self.mdTraj,self.ResidueMasses,array,_cnt_share,lock))
            processes.append(_p)
        
        # start all processes
        for process in processes:
            process.start() 
         
        # wait for all processes to finish
        for process in processes:
            process.join()         
        
        PBar.Clear()
        final = time.time()
        totaltime = final - start
        self.Write2Log('RunTime:            {}\n'.format(totaltime))      
        
        print('\n')
        print('Done Mapping...{0:4.2e} mininutes'.format(totaltime/60.))
        print('Outputting Trajectories...')
        
        # Save the mapped trajectories
        COM_xyz = array
        trajCOM = md.Trajectory(COM_xyz,_top2,unitcell_lengths = self.mdTraj.unitcell_lengths, unitcell_angles = self.mdTraj.unitcell_angles)
        trajCOM.save_lammpstrj(os.path.join(self.SaveDirName,"TrajCOM.lammpstrj"))
        trajCOM.save_dcd(os.path.join(self.SaveDirName,"TrajCOM.dcd"))
        _PDB = md.formats.PDBTrajectoryFile(os.path.join(self.SaveDirName,"TrajCOM.pdb"), mode='w', force_overwrite=True, standard_names=True)
        _PDB.write(trajCOM.xyz[0],trajCOM.topology)
        
        print('Mapping Complete!')