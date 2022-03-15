#!/bin/bash

if [[ -z $1 ]]; then
  echo "ERROR: Command line arguments missing. Syntax: async_pass.sh [input file to be processed] [download mode (optional: default = 0)]"
  exit 1
fi

# check requested files
if [ ! -f $1 ]; then
  echo "ERROR: input file $1 not found locally"
  exit 1
fi
if [ ! -f "o2sim_geometry.root" ]; then
  echo "Cannot find o2sim_geometry.root"
  exit 1
fi
if [ ! -f "o2sim_grp.root" ]; then
  echo "Cannot find o2sim_grp.root"
  exit 1
fi
if [ ! -f "ctf_dictionary.root" ]; then
  echo "Cannot find ctf_dictionary.root"
  exit 1
fi
echo "Checking current directory content:"
ls -altr 

# prepare the list of input files
MODE="local"
if [[ "${1##*.}" == "root" ]]; then
  echo "${1}" > wn.txt
elif [[ "${1##*.}" == "txt" ]]; then
  cp $1 wn.txt
elif [[ "${1##*.}" == "xml" ]]; then
  sed -rn 's/.*turl="([^"]*)".*/\1/p' $1 > wn.txt
  if [ -f `head -n1 wn.txt | xargs -L 1 basename` ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
      sed -i '' -e 's#^.*/##g' wn.txt
    else
      sed -i -e 's#^.*/##g' wn.txt
    fi
  else
    MODE="remote"
  fi
else
  echo "missing input"
  exit 3
fi
echo "processing the following input file(s):"
cat wn.txt

# settings for CTF reader
if [[ $MODE == "remote" ]]; then
  export ARGS_EXTRA_PROCESS_o2_ctf_reader_workflow="--remote-regex \"^alien:///alice/data/.+\""
  if [[ -z $2 || $2 -eq 0 ]]; then
    export INPUT_FILE_COPY_CMD="\"no-copy\""
  else
    export INPUT_FILE_COPY_CMD="\"alien_cp ?src file://?dst\""
  fi
fi

# extra workflow settings
export BEAMTYPE="pp"
export IS_SIMULATED_DATA=0
export DISABLE_MC="--disable-mc"
export DISABLE_ROOT_OUTPUT=""
export TFDELAY=1
export NTIMEFRAMES=-1
export SHMSIZE=16000000000
export WORKFLOW_DETECTORS=MCH,MID
export ARGS_EXTRA_PROCESS_o2_mch_reco_workflow="--digits --triggered"
export CONFIG_EXTRA_PROCESS_o2_mch_reco_workflow="MCHTriggering.triggerRange[0]=-100;MCHTriggering.triggerRange[1]=100;MCHClustering.defaultClusterResolution=0.4;MCHTracking.chamberResolutionX=0.4;MCHTracking.chamberResolutionY=0.4;MCHTracking.sigmaCutForTracking=7.;MCHTracking.sigmaCutForImprovement=6."

# prepare to run workflow
if [[ -f run-workflow-on-inputlist.sh ]]; then
  chmod +x run-workflow-on-inputlist.sh
  echo "Use run-workflow-on-inputlist.sh macro passed as input"
else
  echo "Use run-workflow-on-inputlist.sh macro from O2"
  cp $O2_ROOT/prodtests/full-system-test/run-workflow-on-inputlist.sh .
fi
if [[ -f dpl-workflow.sh ]]; then
  chmod +x dpl-workflow.sh
  echo "Use dpl-workflow.sh macro passed as input"
else
  echo "Use dpl-workflow.sh macro from O2"
  export SETENV_NO_ULIMIT=1
  ln -sf $O2DPG_ROOT/DATA/common/setenv.sh
  cp $O2_ROOT/prodtests/full-system-test/dpl-workflow.sh .
fi

# print workflow
WORKFLOWMODE=print ./run-workflow-on-inputlist.sh CTF wn.txt > workflowconfig.log

# run workflow
WORKFLOWMODE=run ./run-workflow-on-inputlist.sh CTF wn.txt

# check ouput
root -l -b -q 'check.C+({{"mchtracks.root","o2sim"},{"mid-reco.root","midreco"},{"muontracks.root","o2sim"}})'
exit $?
