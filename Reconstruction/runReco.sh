#!/bin/bash
echo ALICE_ROOT = $ALICE_ROOT
echo AliROOT = $AliROOT
cp $ALICE_ROOT/.rootrc ~/.rootrc
cp $ALICE_ROOT/.rootrc $HOME
#cat $HOME/.rootrc
export GRID_TOKEN=OK
export XRD_TRANSACTIONTIMEOUT 300

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

echo ">>>>>>> Running AliRoot to reconstruct $1. Run number is $runnum..."
aliroot -l -b -q rec.C\(\"alien://$1\"\) 2>&1 | tee rec.log
