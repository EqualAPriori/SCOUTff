#import numpy as np
import time
import os
import sys
import re
import subprocess as prcs
import shutil
import parmed
import numpy as np

''' User Inputs '''
nPEO = 10
bPEO = 1  # diblock or triblock
nPPO = 0
nChains = 1 # create two random tacticity chains
tleapName = 'tleap_PEObPPO.in'
itpout = True # split into .itp and .gro file
name = 'NeatPEO_N_10'
PDBName = 'PackmolDir' # file for all the created PDBs for packmol
#Terminate PEO with methyl or hydroxyl
EOmethyl = True

# BUILD PACKMOL:
BuildPackmol = True
NumberConfigs = 120 # total number to pack in
xySideLength = 58. # in Angstroms
zSideLength = 90. # in Angstroms
BoxShell = 0.1 # in Angstroms; keeps molecules away from edge 
Tolerance = 1.25 # in Angstroms
NumberSolvent = 0#10200
SolventPDB = 'water_4.pdb'
FFFile = 'FFFile'


''' CALCULATIONS/FILE SETUP '''
# Create Amber data output directory
saveloc = 'nPEO_{}_nPPO_{}_bPEO_{}_nChains_{}'.format(nPEO,nPPO,bPEO,nChains)
dirName = saveloc
cwd = os.getcwd()

try:
	# Create target Directory
	os.mkdir(os.path.join(cwd,dirName))
	print("Directory {} Created ".format(dirName)) 
except:
	print("WARNING: Overwriting Already Existing Directory {}! ".format(dirName))
	pass

os.mkdir(os.path.join(cwd,dirName,PDBName))
os.mkdir(os.path.join(cwd,dirName,FFFile))
os.chdir(os.path.join(cwd,dirName))

for _chain in range(nChains): # create different chains
	chainName = 'chain_{}'.format(_chain)
	os.mkdir(chainName)
	_name = name +'_{}'.format(_chain) 
	
	print('Chain {}'.format(_chain))
	
	if nPEO == 0: # PPO homopolymer
		nPPO_mid = nPPO - 2
	elif nPEO != 0 and nPPO == 0: # PEO homopolymer 
		nPEO_mid = nPEO - 2
	elif nPEO != 0 and bPEO == 1: # diblock
		nPPO_mid = nPPO - 1
	else: # triblock
		nPPO_mid = nPPO
		
	''' Generate the sequence '''
	if nPPO != 0:
		_seq = np.random.randint(int(0), 2, size=(nPPO_mid), dtype=int)
	
	_seqdict = {0 : 'mP1',
				1 : 'mP2'	
				}
				
	sequence = '{ '
	if nPEO == 0: # homopolymer
		sequence += 'hP1 '
		for _i in range(nPPO-2):
			sequence += '{} '.format(_seqdict[_seq[_i]])
		sequence += 'tP1 }'
	
	elif nPEO != 0 and nPPO == 0: # PEO homopolymer 
		sequence += 'hEO '
		for _i in range(nPEO-1):
			sequence += 'mEO '
		if EOmethyl: # PEO end capping
			sequence += 'tE1 }' # add a methyl
		else:
			sequence += 'tEO }' # add a hydroxyl
	
	elif nPEO != 0 and bPEO == 1: # diblock
		sequence += 'hEO '
		for _i in range(nPEO-1):
			sequence += 'mEO '
		for _i in range(nPPO-1):
			sequence += '{} '.format(_seqdict[_seq[_i]])
		sequence += 'tP1 }'
	
	else: # triblock
		sequence += 'hEO '
		for _i in range(nPEO-1):
			sequence += 'mEO '
		for _i in range(nPPO-2):
			sequence += '{} '.format(_seqdict[_seq[_i]])
		for _i in range(nPEO-1):
			sequence += 'mEO '
		if EOmethyl: # PEO end capping
			sequence += 'tE1 }' # add a methyl
		else:
			sequence += 'tEO }' # add a hydroxyl
	
	print(sequence)
	
	for filename in os.listdir(cwd):
		if filename.endswith(".prepi"):
			shutil.copyfile(os.path.join(cwd,filename), os.path.join(cwd,dirName,chainName,filename))
		else:
			continue
	
	os.chdir(os.path.join(cwd,dirName,chainName))
	
	with open(os.path.join(cwd,tleapName),'r') as tleap:
		ini = tleap.read()
		ini=re.sub('__DUMMYSEQ__',str(sequence),ini)
		ini=re.sub('__NAME__',str(_name),ini)
	
		runfile = open("tleap.in","w")
		runfile.write(ini)
		runfile.close()
	
	call_1 = 'tleap -s -f tleap.in > tleap.out'   

	print(call_1)
	p1 = prcs.Popen(call_1, stdout=prcs.PIPE, shell=True)
	print("p1 \n", p1.communicate())
	
	amber = parmed.load_file(_name+'.parm7', xyz=_name+'.rst7')
	if itpout:
		print("Saving a separate .itp file with the atom,bond types and parameters")
		amber.save(_name+'.top', parameters=_name+'.itp', overwrite=True)
	else:
		amber.save(_name+'.top', overwrite=True)
	amber.save(_name+'.gro', overwrite=True)
	
	# Replace moleculetype/name to chain_#
	with open(_name+'.top','r') as _top:
		ini = _top.read()
		ini=re.sub('system1',str(_name),ini)
	
	os.remove(_name+'.top')
	_top = open(_name+'.top','w')
	_top.write(ini)
	_top.close()
	
	# This is messy, but this reopens the .top file and doesn't write anything else below [ moleucles section ]
	with open(_name+'.top','r') as _f:
		lines = _f.readlines()
	with open(_name+'.top', "w") as _f:
		for line in lines:
			if line.strip("\n") != "[ system ]":
				_f.write(line)
			else: # stop writing anything else out
				break
	
	shutil.copyfile(os.path.join(cwd,dirName,chainName,_name+'.pdb'), os.path.join(cwd,dirName,PDBName,_name+'.pdb'))
	shutil.copyfile(os.path.join(cwd,dirName,chainName,_name+'.top'), os.path.join(cwd,dirName,FFFile,_name+'.top'))
	shutil.copyfile(os.path.join(cwd,dirName,chainName,_name+'.itp'), os.path.join(cwd,dirName,FFFile,_name+'.itp'))
	
	
	os.chdir('..')

''' Build Packmol File '''
if BuildPackmol: 
	os.chdir(PDBName)
	shutil.copyfile(os.path.join(cwd,SolventPDB), os.path.join(cwd,dirName,PDBName,SolventPDB))
	numberperchain = NumberConfigs/nChains
		
	# build packmol 
	Packmol_Input = open('mixture.inp','w')

	Packmol_Input.write('tolerance {}\n'.format(float(Tolerance)))
	Packmol_Input.write('\noutput mixture.pdb\n')
	
	for _chain in range(nChains):
		_name = name +'_{}'.format(_chain) 
		Packmol_Input.write('\nstructure {}\n'.format(_name+'.pdb'))
		Packmol_Input.write('  number {}\n'.format(numberperchain))
		Packmol_Input.write('  inside box {} {} {} {} {} {}\n'.format(BoxShell,BoxShell,BoxShell,(xySideLength-BoxShell),(xySideLength-BoxShell),(zSideLength-BoxShell)))
		Packmol_Input.write('end structure\n')

	# Add Solvent
	if NumberSolvent != 0:
		Packmol_Input.write('\nstructure {}\n'.format(SolventPDB))
		Packmol_Input.write('  number {}\n'.format(NumberSolvent))
		Packmol_Input.write('  inside box {} {} {} {} {} {}\n'.format(BoxShell,BoxShell,BoxShell,(xySideLength-BoxShell),(xySideLength-BoxShell),(zSideLength-BoxShell)))
		Packmol_Input.write('end structure\n')
	
	os.chdir('..')

''' Build the Force-Field file '''
os.chdir(FFFile)
_FF = open('FF.top','w')
_FF.write('; ******* Force-Field File for {} *******\n'.format(dirName))
for _chain in range(nChains):
	_name = name +'_{}'.format(_chain) 
	_FF.write('#include "{}"\n'.format(os.path.join(FFFile,_name+'.itp')))
	_FF.write('#include "{}"\n'.format(os.path.join(FFFile,_name+'.top')))
_FF.write('[ system ]\n')
_FF.write('; Name\n')
_FF.write('Mixture\n')
_FF.write('[ molecules ]\n')
_FF.write('; Compound\n')
for _chain in range(nChains):
	_name = name +'_{}'.format(_chain) 
	_FF.write('{} {}\n'.format(_name,numberperchain))
_FF.close()
	