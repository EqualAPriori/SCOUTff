#!/bin/bash
#Setup a one-species system, with standard parameters
#Example call:
#   ./setup_standard.sh /home/kshen/mylib/packmol 16.0 14pBtrans20.pdb 150 448.15 00

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
TEMP=$5
prefix=$6


mydir=${prefix}.${polymer}_n${npoly}_T${TEMP}
if [ ! -d ${mydir} ]
then
    mkdir $mydir
fi
echo === $mydir ===
cd $mydir

    # Standard parameter values, default
    TEMPANNEAL=773.15

    #frequencies and times below are in terms of picoseconds
    DCDFREQ=100
    THERMOFREQ=10
    ANNEALTIME=1000
    EQTIME=1000
    PRODUCTIONTIME=500000


    #set up script
    sed -e "s/__TEMP__/$TEMP/g" < ${toollib}run_template.py > ./run.py
    sed -i "s/__TEMPANNEAL__/$TEMPANNEAL/g" run.py
    sed -i "s/__DCDFREQ__/$DCDFREQ/g" run.py
    sed -i "s/__THERMOFREQ__/$THERMOFREQ/g" run.py
    sed -i "s/__ANNEALTIME__/$ANNEALTIME/g" run.py
    sed -i "s/__EQTIME__/$EQTIME/g" run.py
    sed -i "s/__PRODUCTIONTIME__/$PRODUCTIONTIME/g" run.py


    #auto-detect forcefield // make sure forcefields are not redundant (e.g. re-defining atom types)

    #suggested commands
    sed -e "s/__label__/$mydir/" < ${toollib}podsubmit_template.sh > ./podsubmit.sh

    echo "# Suggested commands" > z.commands
    echo "python ${toollib}make_topology.py -L $L -pkml $pkmldir -structlib $structlib -m $polymer $npoly -ff ${fflib}TrappeUA_Master.xml" -xyz >> z.commands
    echo "srun --gres=gpu:1" "python ${toollib}make_topology.py -L $L -pkml $pkmldir -structlib $structlib -m $polymer $npoly -ff ${fflib}TrappeUA_Master.xml" -xyz >> podsubmit.sh
    echo "" >> z.commands

    echo "python ./run.py -ff ${fflib}TrappeUA_Master.xml -PME -tail" >> z.commands
    echo "srun --gres=gpu:1" "python ./run.py -ff ${fflib}TrappeUA_Master.xml -PME -tail" >> podsubmit.sh
    echo "" >> z.commands
    echo "#python ./run.py -prefix run1 -sysxml run0_system.xml -chkstate run0_checkpoint.xml -PME -tail" >> z.commands
    echo "#srun --gres=gpu:1 python ./run.py -prefix run1 -sysxml run0_system.xml -chkstate run0_checkpoint.xml -PME -tail" >> podsubmit.sh

cd ..
