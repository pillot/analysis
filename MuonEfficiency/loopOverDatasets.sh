#!/bin/sh

# usage:
# $1 = file containing the datasets


cat $1 | while read dataset; do
  
  run=${dataset:(-6)};
  echo "processing run $run..."
  
  mkdir -p runs/$run
  cd runs/$run
  
  root -q -l /Users/pillot/Work/Alice/Work/Macros/MuonEfficiency/runMuonEfficiency.C\(\""saf"\",\"$dataset\"\) >& run.log
  
  cd ../..
  
done

