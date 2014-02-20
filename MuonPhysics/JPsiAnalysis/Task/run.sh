#!/bin/bash

source /Users/lardeux/.bash_profile

export Study=SL
mkdir ./Outputs/$Study

# Task Mixing
export Filename=ListDataSet_132runs_AOD101.txt

root -l << EOF
gEnv->SetValue("XSec.GSI.DelegProxy", "2")
.x RunALICE.C("proof","full","aod","",1e10,0,"AddAMEventMixingTest","$ALICE_ROOT","",kTRUE,"$Filename")
cin>>k
.q
EOF


echo ""
echo "######################################################"
echo "#                 Script is done !                   #"
echo "######################################################"