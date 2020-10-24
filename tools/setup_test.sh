#!/bin/bash
#Setup a one-species system, with standard parameters
#Example call:
#   ./setup_test.sh /home/kshen/mylib/packmol 10.0 14pBtrans20.pdb 10

# Get directory where SCOUTff is and some shortcuts
thisScriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
root=${thisScriptDir}/../
structlib=${root}structlib/Kraton/
toollib=${root}tools/
fflib=${root}TrappeUA/XML/


pkmldir=$1
L=$2
polymer=$3
npoly=$4
TEMP=448.15
prefix=quicktest


mydir=${prefix}.${polymer}_n${npoly}_T${TEMP}
if [ ! -d ${mydir} ]
then
    mkdir $mydir
fi
echo === $mydir ===
cd $mydir

    # Standard parameter values, default
    TEMPANNEAL=558.15

    #frequencies and times below are in terms of picoseconds
    DCDFREQ=10
    THERMOFREQ=1
    ANNEALTIME=10
    EQTIME=10
    PRODUCTIONTIME=10


    #set up script
    sed -e "s/__TEMP__/$TEMP/g" < ${toollib}run_template.py > ./run.py
    sed -i "s/__TEMPANNEAL__/$TEMPANNEAL/g" run.py
    sed -i "s/__DCDFREQ__/$DCDFREQ/g" run.py
    sed -i "s/__THERMOFREQ__/$THERMOFREQ/g" run.py
    sed -i "s/__ANNEALTIME__/$ANNEALTIME/g" run.py
    sed -i "s/__EQTIME__/$EQTIME/g" run.py
    sed -i "s/__PRODUCTIONTIME__/$PRODUCTIONTIME/g" run.py
    #sed -i "0,/use_gpu.*/s/\(use_gpu\).*/\1 = false/" run.py


    #auto-detect forcefield // make sure forcefields are not redundant (e.g. re-defining atom types)

    #suggested commands
    echo "# Suggested commands" > z.commands
    echo "python ${toollib}make_topology.py -L $L -pkml $pkmldir -structlib $structlib -m $polymer $npoly -ff ${fflib}TrappeUA_Master.xml" >> z.commands
    echo "" >> z.commands
    echo "python ./run.py -ff ${fflib}TrappeUA_Master.xml" >> z.commands

cd ..
