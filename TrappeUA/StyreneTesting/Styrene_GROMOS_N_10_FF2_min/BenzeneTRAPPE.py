from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm import app
from sys import stdout
import mdtraj as md
import numpy as np
import parmed
import os
import sys
import re
import subprocess as prcs
import shutil
import time

#First we define utility functions to get information about a system
def groupForces(sysObj):
	#Assign each for a name and an index
	#Returns a dictionary with force names as keys and force indices in the system object as definitions
	# i.e. HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce, NonbondedForce, CMMotionRemover
	forces = {}
	forcesCounts = {}
	for k, frc in enumerate(sysObj.getForces()):
		frcName = frc.__class__.__name__
		if frcName in forces.keys():
			forcesCounts[frcName] += 1
			frcName += '_%i'%forcesCounts[frcName]
		else:
			forcesCounts[frcName] = 1
		forces[frcName] = k
	return forces

logname = 'log.out'
def write2log(logname,text):
    try:
        f = open(logname,'a')
    except:
        f = open(logname,'w')
    f.write(text)
    f.close()

''' Simulation Settings '''
Pressure = 1.0 # in units of bar
Temperature = 463 # units kelvin
barostatfreq = 20
BoxL = 50.0 # Angstrom

NPol = 40
DOP = 10

IncludeHydrogens = False
Toulene = True
BuildPackmol = True
RunNPT = True
TestBenzeneGeom = True
AddContraints = False
useGPU = True
PS = True

''' --- the actual work of creating the topology --- '''
top = app.topology.Topology()

chains = []
residues = []
constrainedpairs = []
atomsbackbone = []

if PS == False: 
    # add Benzene
    for _i in range(NPol): # quick-n-dirty Benzene topology
        chain = top.addChain() #Create new chain for each molecule
        chains.append(chain)
        resname = 'BEN'
        res = top.addResidue(resname,chains[-1])
        residues.append(res)
        
        _C = app.element.Element.getBySymbol("C")
        _H = app.element.Element.getBySymbol("H")
        a1 = top.addAtom('C1', _C, residues[-1] )
        a2 = top.addAtom('C2', _C, residues[-1] )
        a3 = top.addAtom('C3', _C, residues[-1] )
        a4 = top.addAtom('C4', _C, residues[-1] )
        a5 = top.addAtom('C5', _C, residues[-1] )
        a6 = top.addAtom('C6', _C, residues[-1] )
        
        top.addBond( a1, a2 )
        top.addBond( a2, a3 )
        top.addBond( a3, a4 )
        top.addBond( a4, a5 )
        top.addBond( a5, a6 )
        top.addBond( a6, a1 )
        
        if Toulene: # add toulene
            a7 = top.addAtom('C7', _C, residues[-1] )
            top.addBond( a4, a7 )
        
        _temppairs = [[a1.index,a2.index],
                      [a2.index,a3.index],
                      [a3.index,a4.index],
                      [a4.index,a5.index],
                      [a5.index,a6.index],
                      [a6.index,a1.index],
                      [a4.index,a7.index]]
                      
        constrainedpairs.extend(_temppairs) 
        
        if IncludeHydrogens:
            a7 = top.addAtom('H1', _H, residues[-1] )
            a8 = top.addAtom('H2', _H, residues[-1] )
            a9 = top.addAtom('H3', _H, residues[-1] )
            a10 = top.addAtom('H4', _H, residues[-1] )
            a11 = top.addAtom('H5', _H, residues[-1] )
            a12 = top.addAtom('H6', _H, residues[-1] )
            if Toulene:
                a13 = top.addAtom('H7', _H, residues[-1] )
                a14 = top.addAtom('H8', _H, residues[-1] )
                a15 = top.addAtom('H9', _H, residues[-1] )
            
            top.addBond( a1, a7 )
            top.addBond( a2, a8 )
            top.addBond( a3, a9 )
            top.addBond( a4, a10 )
            top.addBond( a5, a11 )
            top.addBond( a6, a12 )
            if Toulene:
                top.addBond( a4, a13 )
                top.addBond( a4, a14 )
                top.addBond( a4, a15 )
elif PS:
    # add Polystrene
    for _i in range(NPol): # quick-n-dirty PS topology
        chain = top.addChain() #Create new chain for each molecule
        chains.append(chain)
        for _j in range(DOP):
            resname = 'PS'
            res = top.addResidue(resname,chains[-1])
            residues.append(res)
            
            # add the atoms
            if _j == 0: # first monomer
                _C = app.element.Element.getBySymbol("C")
                _O = app.element.Element.getBySymbol("O")
                _H = app.element.Element.getBySymbol("H")
                a1 = top.addAtom('C1', _C, residues[-1] )
                a2 = top.addAtom('C2', _C, residues[-1] )
                a3 = top.addAtom('C3', _C, residues[-1] )
                a4 = top.addAtom('C4', _C, residues[-1] )
                a5 = top.addAtom('C5', _C, residues[-1] )
                a6 = top.addAtom('C6', _C, residues[-1] )
                a7 = top.addAtom('C7', _C, residues[-1] )
                a8 = top.addAtom('C8', _C, residues[-1] )
 
                top.addBond( a1, a2 )
                top.addBond( a2, a3 )
                top.addBond( a3, a4 )
                top.addBond( a4, a5 )
                top.addBond( a5, a6 )
                top.addBond( a6, a1 )
                top.addBond( a7, a1 )
                top.addBond( a7, a8 )            

                atomsbackbone.append([a7,a8])
            
            elif _j == int(DOP-1): # last monomer
                _C = app.element.Element.getBySymbol("C")
                _O = app.element.Element.getBySymbol("O")
                _H = app.element.Element.getBySymbol("H")
                a1 = top.addAtom('C1', _C, residues[-1] )
                a2 = top.addAtom('C2', _C, residues[-1] )
                a3 = top.addAtom('C3', _C, residues[-1] )
                a4 = top.addAtom('C4', _C, residues[-1] )
                a5 = top.addAtom('C5', _C, residues[-1] )
                a6 = top.addAtom('C6', _C, residues[-1] )
                a7 = top.addAtom('C7', _C, residues[-1] )
                a8 = top.addAtom('C8', _C, residues[-1] )
                
                top.addBond( a1, a2 )
                top.addBond( a2, a3 )
                top.addBond( a3, a4 )
                top.addBond( a4, a5 )
                top.addBond( a5, a6 )
                top.addBond( a6, a1 )
                top.addBond( a7, a1 )
                top.addBond( a7, a8 )              
                   
                top.addBond( a7, atomsbackbone[-1][-1])
                atomsbackbone.append([a7,a8]) 
            
            else: # middle monomers
                _C = app.element.Element.getBySymbol("C")
                _O = app.element.Element.getBySymbol("O")
                _H = app.element.Element.getBySymbol("H")

                a1 = top.addAtom('C1', _C, residues[-1] )
                a2 = top.addAtom('C2', _C, residues[-1] )
                a3 = top.addAtom('C3', _C, residues[-1] )
                a4 = top.addAtom('C4', _C, residues[-1] )
                a5 = top.addAtom('C5', _C, residues[-1] )
                a6 = top.addAtom('C6', _C, residues[-1] )
                a7 = top.addAtom('C7', _C, residues[-1] )
                a8 = top.addAtom('C8', _C, residues[-1] )
                
                top.addBond( a1, a2 )
                top.addBond( a2, a3 )
                top.addBond( a3, a4 )
                top.addBond( a4, a5 )
                top.addBond( a5, a6 )
                top.addBond( a6, a1 )
                top.addBond( a7, a1 )
                top.addBond( a7, a8 )  

                top.addBond( a7, atomsbackbone[-1][-1])
                atomsbackbone.append([a7,a8])

element_VS = app.element.Element(201, "M", "VS", 0.*dalton)    

# create force field 
# load in forcefield
forcefield = ForceField('TrappeUA_Styrene_Gromos.xml','tip4pew.xml')
# look for unmatched residues
unmatched_residues = forcefield.getUnmatchedResidues(top)
write2log(logname,"Unmatched\n")
if len(unmatched_residues) > 0:
    write2log(logname,"{} \n".format(unmatched_residues))
else:
    write2log(logname,"{} \n".format('None'))

# Get positions
if BuildPackmol == False:
    Pos = np.random.random((top.getNumAtoms(),3))*BoxL       
    print("Position Shape: {}".format(np.shape(Pos)))

periodicboxvectors = [[BoxL/10.,0.,0.],[0.,BoxL/10.,0.],[0.,0.,BoxL/10.]]
top.setPeriodicBoxVectors(periodicboxvectors*nanometer)
   
''' Build Packmol File '''
if BuildPackmol: 
    cwd = os.getcwd()
    PackmolDir = 'packmol'
    PDB = 'loadpdb_np_01_N_10.pdb'
    BoxShell = 0.1
    Tolerance = 2.5
    try:
        os.mkdir(PackmolDir)
    except:
        shutil.rmtree(PackmolDir)
        os.mkdir(PackmolDir)
        
    os.chdir(PackmolDir)
    shutil.copyfile(os.path.join(cwd,PDB), os.path.join(cwd,PackmolDir,PDB))
        
    # build packmol 
    Packmol_Input = open('mixture.inp','w')

    Packmol_Input.write('tolerance {}\n'.format(float(Tolerance)))
    Packmol_Input.write('\noutput mixture.pdb\n')

    Packmol_Input.write('\nstructure {}\n'.format(PDB))
    Packmol_Input.write('  number {}\n'.format(NPol))
    Packmol_Input.write('  inside box {} {} {} {} {} {}\n'.format(BoxShell,BoxShell,BoxShell,(BoxL*10-BoxShell),(BoxL*10-BoxShell),(BoxL*10-BoxShell)))
    Packmol_Input.write('end structure\n')
    Packmol_Input.close()

    call_1 = 'packmol < mixture.inp'   

    print(call_1)
    p1 = prcs.Popen(call_1, stdout=prcs.PIPE, shell=True)
    print("p1:\n {}".format(p1.communicate()))

    shutil.copyfile(os.path.join(cwd,PackmolDir,'mixture.pdb'), os.path.join(cwd,'mixture.pdb'))
    os.chdir('..')
   
''' Create System '''
system = forcefield.createSystem(top, 
                                nonbondedMethod=LJPME,
                                nonbondedCutoff=1.*nanometer, 
                                ewaldErrorTolerance=0.0001, 
                                rigidWater=True, 
                                constraints=None)
                                
integrator = LangevinIntegrator(Temperature*kelvin, 1/picosecond, 0.002*picoseconds)
write2log(logname,'Integrator Constraint Tolerance: {}'.format(integrator.getConstraintTolerance()))
barostat = MonteCarloBarostat(Pressure*bar, Temperature*kelvin, barostatfreq)
if RunNPT: # add barostat
    write2log(logname,"Adding Barostat \n")
    system.addForce(barostat)

if useGPU:
	platform = Platform.getPlatformByName('OpenCL') 
	properties = {'DeviceIndex': '1', 'Precision': 'mixed'}
else:
	platform = Platform.getPlatformByName('CPU')
	properties = {'Threads': '10'}

# add benzene ring constraints
if AddContraints:
    write2log(logname,"Shape of constrained pairs: {} \n".format(np.shape(constrainedpairs)))
    for _i,pair in enumerate(constrainedpairs):
        length = 0.140
        system.addConstraint(pair[0],pair[1],length*nanometer)

# setup simulation
simulation = Simulation(top, system, integrator, platform, properties)
if BuildPackmol:
    pdb = PDBFile('mixture.pdb')
    positions = pdb.positions # these are in nanometers
else:
    positions = Pos*nanometer

simulation.context.setPositions(positions)
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('test.pdb','w'))
simulation.context.setPeriodicBoxVectors(periodicboxvectors[0],periodicboxvectors[1],periodicboxvectors[2]) # Set the periodic box vectors

# write out xml
from simtk.openmm import XmlSerializer
serialized_system_gromacs = XmlSerializer.serialize(system)
outfile = open('system.xml','w')
outfile.write(serialized_system_gromacs)
outfile.close()

# report on platform
periodicboxvectors = simulation.context.getState().getPeriodicBoxVectors() # Get the periodic box vectors	
p = simulation.context.getPlatform()
force_dict = groupForces(system)
nbf_obj = system.getForces()[force_dict['NonbondedForce']]

write2log(logname,"Box Size: {}".format(periodicboxvectors))
write2log(logname,"simulation platform: {} \n".format(p.getName()) )
write2log(logname,"{} \n".format(p.getPropertyNames()))
if useGPU:
    write2log(logname,"Device name: {} \n".format(p.getPropertyValue(simulation.context,'DeviceName')))
    write2log(logname,"Use CPU PME: {} \n".format(p.getPropertyValue(simulation.context,'UseCpuPme')))
    write2log(logname,"Disable PME Stream: {} \n".format(p.getPropertyValue(simulation.context,'DisablePmeStream')))
    write2log(logname,"Device Index: {} \n".format(p.getPropertyValue(simulation.context,'DeviceIndex')))
write2log(logname,"Nonbonded method: {} \n".format(nbf_obj.getNonbondedMethod()))
write2log(logname,'The dispersion correction setting is: {}.\n'.format(nbf_obj.getUseDispersionCorrection()))

''' Save bond parameters to file. '''
numBonds = system.getForces()[force_dict['HarmonicBondForce']].getNumBonds()
with open("bond.parameters", 'w') as g:
	for i in range(numBonds):
		g.write('bond {}: {} \n'.format(i,system.getForces()[force_dict['HarmonicBondForce']].getBondParameters(i)))
		
''' Save angle parameters to file. '''
numAngles = system.getForces()[force_dict['HarmonicAngleForce']].getNumAngles()
with open("Angle.parameters", 'w') as g:
	for i in range(numAngles):
		g.write('angle {}: {} \n'.format(i,system.getForces()[force_dict['HarmonicAngleForce']].getAngleParameters(i)))
	
''' Save torsion parameters to file. '''
numTorsions = system.getForces()[force_dict['CustomTorsionForce']].getNumTorsions()
with open("Torsion.parameters", 'w') as g:
	for i in range(numTorsions):
		g.write('torsion {}: {} \n'.format(i,system.getForces()[force_dict['CustomTorsionForce']].getTorsionParameters(i)))
		
''' Save non-bonded parameters to file. '''
numAtoms = nbf_obj.getNumParticles()
with open("nonbonded.parameters", "w") as g:
	for i in range(numAtoms):
		param = nbf_obj.getParticleParameters(i)
		''' Modify non-bonded forces if desired. '''
		#nbf_obj.setParticleParameters(index=i, charge=scalecharge*param[0].value_in_unit(elementary_charge), sigma=param[1].value_in_unit(nanometer), epsilon=param[2].value_in_unit(kilojoules_per_mole) )
		#param_updated = nbf_obj.getParticleParameters(i)
		mass = system.getParticleMass(i)
		g.write('atom {} and mass {}: {}\n'.format(i,mass,param))

# minimize and run simulation
simulation.context.applyConstraints(1e-8)
positions = simulation.context.getState(getPositions=True).getPositions()


# get forces on each particle
forces = simulation.context.getState(getForces=True).getForces()
flog = open('forces.out', 'w')
for i, f in enumerate(forces):
    if norm(f) > 1e4*kilojoules_per_mole/nanometer:
        flog.write("{} {}\n".format(i,f))
flog.close()

write2log(logname,"Done with Constraints \n")
timestart = time.time()
write2log(logname,'pre minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy(tolerance=0.01*kilojoules_per_mole,maxIterations=1000000)
write2log(logname,'post minimization potential energy: {} \n'.format(simulation.context.getState(getEnergy=True).getPotentialEnergy()))
timeend = time.time()
write2log(logname,"Done with Minimization in {} mininutes\n".format((timeend-timestart)/60.))
simulation.reporters.append(app.dcdreporter.DCDReporter('npt.dcd', 1000))
simulation.reporters.append(StateDataReporter('thermo.out', 1000, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True))

forces = simulation.context.getState(getForces=True).getForces()
flog = open('forces_1.out', 'w')
for i, f in enumerate(forces):
    if norm(f) > 1e4*kilojoules_per_mole/nanometer:
        flog.write("{} {}\n".format(i,f))
flog.close()

simulation.step(5000) # run short NVT
write2log(logname,"Done with Short NVT \n")

if TestBenzeneGeom: # check if ring planar
    gout = open('geom.out','w')
    gout.write('#  d25  d36  d47  |  d37  ||  d25-d36  d25-d47  d36-d47  |||  chk1  chk2  control\n')
    gout.close()
    for _i in range(100):
        positions = simulation.context.getState(getPositions=True).getPositions()
        d25 = norm(positions[0]-positions[3])
        d36 = norm(positions[1]-positions[4])
        d47 = norm(positions[2]-positions[5])
        d37 = norm(positions[1]-positions[5])
        pvec = np.cross(positions[0],positions[3])
        chk1 = np.inner(positions[2],pvec)
        chk2 = np.inner(positions[4],pvec)
        chk3 = np.inner(positions[0],pvec)
        gout = open('geom.out','a')
        gout.write('{}  {}  {}  {}  |  {}  ||  {}  {}  {}  |||  {}  {}  {}\n'.format(_i,d25,d36,d47,d37,d25-d36,d25-d47,d36-d47,chk1,chk2,chk3))
        gout.close()
        simulation.step(200)
        

positions = simulation.context.getState(getPositions=True).getPositions()       
app.pdbfile.PDBFile.writeModel(simulation.topology,positions,open('loadpdb.pdb','w'))
simulation.step(50E6)