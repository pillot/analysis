#!/bin/sh

# usage:
# $1 = mode
# $2 = dataset

if [ $# != 2 ]; then
    echo "usage: \$1=\"local\", \"saf\" or \"saf3\"; \$2=dataset"
    exit 1
fi

runonsaf3="no";

if [ "$1" = "saf3" ] && [ `hostname` = "nansafmaster3.in2p3.fr" ]; then

    vafctl start
    nWorkers=88
    let "nWorkers -= `pod-info -n`"
    echo "requesting $nWorkers additional workers"
    vafreq $nWorkers
    vafwait 88
    runonsaf3="yes"

fi

for (( ipT=0 ; ipT < 4 ; ipT++ )); do

    for (( iy=0 ; iy < 12 ; iy++ )); do
	
        fileName="acceff_pT${ipT}_y${iy}.root"

        echo "processing pT$ipT - y$iy..."

        if [ -e $fileName ]; then
            continue
        fi

        if [ "$runonsaf3" = "yes" ]; then
            root -q -b runJPsiAccEffCorr2.C\(\"$1\",\"$2\",$ipT,$iy\) 2>&1 | tee run_pT${ipT}_y${iy}.log
        else
            root -q -l $WORK/Macros/MuonEfficiency/runJPsiAccEffCorr2.C\(\"$1\",\"$2\",$ipT,$iy\) 2>&1 | tee run_pT${ipT}_y${iy}.log
        fi

        mv acceff_new.root $fileName

    done

done

if [ "$runonsaf3" = "yes" ]; then

    vafctl stop

fi
