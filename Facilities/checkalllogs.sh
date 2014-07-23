#!/bin/sh

#usage:
#$1 = runList
#$2 = outDir

if [[ ! -e $1 ]]; then
  echo "provide the run list"
  exit 2
fi

outdir="CheckLogs"
if [[ $2 ]]; then
  outdir=$2
fi

if [[ -d $outdir ]]; then
  echo "a directory containing checks ($outdir/) already exist. Please move or remove it"
  exit 3
fi

mkdir $outdir

for run in `cat $1`; do
  
  cd $outdir
#  $WORK/Macros/Facilities/checklogs.sh /alice/data/2011/LHC11h/000$run/ESDs/pass2_muon $run &
  $WORK/Macros/Facilities/checklogs.sh /alice/data/2012/LHC12h/000$run/ESDs/muon_calo_pass2 $run &
  cd ..
  
done
