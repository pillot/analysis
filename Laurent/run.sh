#!/bin/bash

#source /Users/pillot/Work/bin/alienenv.sh /Users/pillot/Work/Alice/Work

rm -f results.list

for ((i=1; i <= 10 ; i++))
do
  root -l << EOF
  .x $ALICE/Macros/Laurent/run.C("saf", "../AOD49.shortlist$i", "AnalysisResults_part$i.root", 'a')
  .q
EOF
  echo "AnalysisResults_part$i.root" >> results.list
done

root -l -b << EOF
gROOT->LoadMacro("$ALICE/Macros/Facilities/runTaskFacilities.C");
LoadAlirootLocally("PWG3base", "", "AliHistogramCollection");
.x $ALICE_ROOT/PWG3/muon/mergeGridFiles.C("AnalysisResults.root", "results.list", "", 1000, kFALSE, "");
.q
EOF
