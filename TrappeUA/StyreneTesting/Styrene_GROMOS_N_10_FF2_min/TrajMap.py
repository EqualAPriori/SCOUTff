'''
Example Script for Using the TrajMap_Module.

Can use both serial and pll to process a trajectory. Currently only setup to map to COM of residues in PDB.

PLL Usage:
    - use PLL with ~300-500 residues per thread or more. Otherwise serial is probably fine. 

'''

import TrajMap_Module as TM

warmup = 250
stride = 1

DoResidueBasedMapping = False
DoResidueBasedMappingPLL = True # uses threading
_nProcessors = 6
# ToDo
DoAtomBasedMapping = False
DoCustomResidueMapping = False

# Does COM mapping based off of residues in mixture.pdb
if DoResidueBasedMapping:
    TM = TM.TrajObj(['npt.dcd'],'initial.pdb',
                    _Warmup=warmup,_Stride=stride,_CreateSaveDir=True,_SaveDirName='TrajMap_ResidueBased_Wrapped')
    TM.DoCOMMap()

# Does COM mapping based off of residues in mixture.pdb but uses threading
if DoResidueBasedMappingPLL:
    TM = TM.TrajObj(['npt.dcd'],'loadpdb.pdb',
                    _Warmup=warmup,_Stride=stride,_CreateSaveDir=True,_SaveDirName='TrajMap_ResidueBasedPLL_np_{}'.format(_nProcessors))
    TM.DoCOMMapPLL(_nProcessors=_nProcessors)

# TODO: Not implemented yet
if DoAtomBasedMapping:
    TM = TM.TrajObj(['nvt_melt.dcd'],'mixture.pdb',
                    _Warmup=warmup,_Stride=stride,_CreateSaveDir=True,_SaveDirName='TrajMap_AtomBased')
    TM.DoCOMMap()
    
if DoCustomResidueMapping:
    TM = TM.TrajObj(['nvt_melt.dcd'],'mixture.pdb',
                    _Warmup=warmup,_Stride=stride,_CreateSaveDir=True,_SaveDirName='TrajMap_CustomResidue')
    TM.DoCOMMap()