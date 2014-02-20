#!/bin/sh

# usage:
# $1 = dataset
# $2 = tracker or trigger

if [ "$2" = "tracker" ]; then
    
    sigmaTrg=4
    
#    nTest=7
#    resoNB=(0.75 0.75 1 1 1 2 4)
#    resoB=(0.25 0.5 0.5 0.75 1 2 4)
    nTest=8
    resoNB=(0.75 0.75 1 1 1.5 2 3 4)
    resoB=(0.25 0.5 0.75 1 1.5 2 3 4)
    
    for (( i=0 ; i < nTest ; i++ )); do
	
	errB=`echo ${resoB[i]} | sed 's/0\./0\.0/g'`
	errB=`echo ${errB} | sed 's/1\./0\.1/g'`
	if [ "$errB" = "${resoB[i]}" ]; then errB="0.${resoB[i]}"; fi
	errNB=`echo ${resoNB[i]} | sed 's/0\./0\.0/g'`
	errNB=`echo ${errNB} | sed 's/1\./0\.1/g'`
	if [ "$errNB" = "${resoNB[i]}" ]; then errNB="0.${resoNB[i]}"; fi
	
	for sigma in 2 3 4 5 6 ; do
#	for sigma in 2 3 4 ; do
	    
	    dirname="${resoNB[i]}-${resoB[i]}mm_${sigma}sigma_${sigmaTrg}sigma"
	    
	    echo "processing $dirname..."
	    
	    if [ -e $dirname ]; then
		continue
	    fi
	    
	    mkdir $dirname
	    
	    cd $dirname
	    
	    root -q -l /Users/pillot/Work/Alice/Work/Macros/Refit/runMuonRefit.C\(\""saf"\",\"$1\",$errNB,$errB,$sigma,$sigmaTrg\) >& run.log
	    
	    cd ..
	    
	done
	
    done
    
elif [ "$2" = "trigger" ]; then
    
    resoNB=2
    resoB=2
    sigmaTrk=4
    
    errB=`echo ${resoB[i]} | sed 's/\./\.0/g'`
    if [ "$errB" = "${resoB[i]}" ]; then errB="0.${resoB[i]}"; fi
    errNB=`echo ${resoNB[i]} | sed 's/\./\.0/g'`
    if [ "$errNB" = "${resoNB[i]}" ]; then errNB="0.${resoNB[i]}"; fi
    
    for sigma in 2 3 4 5 6 ; do
#    for sigma in 2 3 4 ; do
#    for sigma in 4 ; do
	
	dirname="${resoNB}-${resoB}mm_${sigmaTrk}sigma_${sigma}sigma"
	
	echo "processing $dirname..."
	
	if [ -e $dirname ]; then
	    continue
	fi
	
	mkdir $dirname
	
	cd $dirname
	
	root -q -l /Users/pillot/Work/Alice/Work/Macros/Refit/runMuonRefit.C\(\""saf"\",\"$1\",$errNB,$errB,$sigmaTrk,$sigma\) >& run.log
#	root -q -l /Users/pillot/Work/Alice/Work/Macros/Refit/runMuonRefit.C\(\""local"\",\"$1\",$errNB,$errB,$sigmaTrk,$sigma\) >& run.log
	
	cd ..
	
    done
    
else
    
    echo "usage: \$1=dataset; \$2=\"tracker\" or \"trigger\""
    
fi

