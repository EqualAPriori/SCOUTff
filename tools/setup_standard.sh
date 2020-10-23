#!/bin/bash
#Setup a one-species system, with standard parameters
#Example call:
#   ../setup_standard.sh 10.0 14pBtrans20.pdb 10 ~/openmm/Structures/packmol

L=$1
polymer=$2
npoly=$3
pkmldir=$4

# Get directory where SCOUTff is and some shortcuts
thisScriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
root=${thisScriptDir}/../
structlib=${root}structlib/Kraton/
toollib=${root}tools/
fflib=${root}TrappeUA/XML/

# Standard parameter values, default
TEMP=448.15
TEMPANNEAL=558.15

#frequencies and times below are in terms of picoseconds
DCDFREQ=100
THERMOFREQ=1
ANNEALTIME=1000
EQTIME=1000
PRODUCTIONTIME=100000


#set up script
sed -e "s/__TEMP__/$TEMP/g" < ${toollib}/run_template.py > ./run.py
sed -i "s/__TEMPANNEAL__/$TEMPANNEAL/g" run.py
sed -i "s/__DCDFREQ__/$DCDFREQ/g" run.py
sed -i "s/__THERMOFREQ__/$THERMOFREQ/g" run.py
sed -i "s/__ANNEALTIME__/$ANNEALTIME/g" run.py
sed -i "s/__EQTIME__/$EQTIME/g" run.py
sed -i "s/__PRODUCTIONTIME__/$PRODUCTIONTIME/g" run.py


#auto-detect forcefield // make sure forcefields are not redundant (e.g. re-defining atom types)

#suggested commands
echo "# Suggested commands" > z.commands
echo "python ${toollib}make_topology.py -L $L -pkml $pkmldir -structlib $structlib -m $polymer $npoly -ff ${fflib}TrappeUA_Butadiene_Gromos.xml" >> z.commands
echo "" >> z.commands
echo "python ./run.py -ff ${fflib}TrappeUA_Butadiene_Gromos.xml" >> z.commands


