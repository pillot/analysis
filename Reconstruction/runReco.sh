#!/bin/bash
echo ALICE_ROOT = $ALICE_ROOT
echo AliROOT = $AliROOT
cp $ALICE_ROOT/.rootrc ~/.rootrc
cp $ALICE_ROOT/.rootrc $HOME
cat $HOME/.rootrc
export GRID_TOKEN=OK

echo ">>>>>>>>> PATH is..."
echo $PATH
echo ">>>>>>>>> LD_LIBRARY_PATH is..."
echo $LD_LIBRARY_PATH
echo ">>>>>>>>> rec.C is..."
cat rec.C
echo

ls -l

# $1 = raw input filename
runnum=`echo $1 | cut -d "/" -f 6`

echo ">>>>>>> Running AliRoot to reconstruct $1 with option \"$2\". Run number is $runnum..."
aliroot -l -b -q runReco.C\(\""alien://$1$2"\"\) 2>&1 | tee rec.log

ls -l | tee -a rec.log
