#!/bin/bash

scriptName=`basename $0`

################## Please change me! #####################
#subDir="outer";
subDir="";
##########################################################

qaFileName=${1:-"QAresults_outer.root"};
isBarrel=`echo "$qaFileName" | grep -c "barrel"`
if [ $isBarrel -gt 0 ]; then
    subDir="barrel"
fi

##########################
# Load aliroot environment
##########################
if [ ! $ALICE_ROOT ]; then
    if [ ! $WORK ]; then
        source $HOME/.profile
    fi
    source /Users/pillot/Work/bin/alienenv.sh /Users/pillot/Work/Alice/trunk2
fi

################## Please change me! #####################
qaPass="pass1" #"pass1_muon" #"pass1_HLT" #pass1_HLT/QA90
#qaPass="pass2_muon" #"pass1_muon" #"pass1_HLT" #pass1_HLT/QA90
alienBaseDir="alien:///alice/data/2012/LHC12b"
#alienBaseDir="alien:///alice/data/2011/LHC11h"
outDir="/Users/pillot/Work/Alice/Work/Data/2012/LHC12b/MuonQA/$qaPass/$subDir"
#outDir="/Users/pillot/Work/Alice/Work/Data/2011/LHC11h/$qaPass/MuonQA/$subDir"
userMail="pillot@subatech.in2p3.fr" #"stocco@subatech.in2p3.fr" #" cynthia@ipno.in2p3.fr palash.khan@cern.ch"

#prodUrlInMonalisa="http://alimonitor.cern.ch/prod/jobs.jsp?t=2656" # 2304 1819" 1830 # If url is set, get runList from url
#prodUrlInMonalisa="http://alimonitor.cern.ch/prod/jobs.jsp?t=1819" #LHC11h
prodUrlInMonalisa="http://alimonitor.cern.ch/prod/jobs.jsp?t=2388" # LHC12a/b
#prodUrlInMonalisa="http://alimonitor.cern.ch/prod/jobs.jsp?t=2698" # LHC12e
#prodUrlInMonalisa="http://alimonitor.cern.ch/prod/jobs.jsp?t=2715" # LHC12f
alirootInMonalisa="v5-02-Rev-16a" #"v5-01-Rev" #"v5-01-Rev-12" #"v5-02-09-AN"
qaPassInMonalisa="$qaPass"

alienUserName="${alien_API_USER}"
alienCloseSE="${alien_CLOSE_SE}"

outAlienDir="alien:///alice/cern.ch/user/${alienUserName:0:1}/${alienUserName}/Data/LHC12b/$qaPass/QA"
#outAlienDir="alien:///alice/cern.ch/user/${alienUserName:0:1}/${alienUserName}/Data/LHC12f/QA/$subDir"
#outAlienDir="alien:///alice/cern.ch/user/${alienUserName:0:1}/${alienUserName}/Data/LHC11h/$qaPass/QA"
runListName="runList.txt"
triggerListName="triggerList.txt"
excludeRunListName="excludeRunList.txt"

outFilename="$outDir/outCronQA.txt"
##########################################################

###########################################
# Perform checks each time asked by crontab
# but re-launch in case of problems only.
# If previous run was successfull
# relaunch only after "waitingTime" hours
###########################################
checkTime=`date "+%y%m%d%H"`
waitingTime=6 # hours from successfull submission
if [ -e $outFilename ]; then
    oldCheckTime=`grep CheckTime ${outFilename} | cut -d " " -f 2`
    if [ "${oldCheckTime}" ]; then
	let "timeDiff = ${checkTime} - ${oldCheckTime}"
	if [ ${timeDiff} -lt ${waitingTime} ]; then
	    echo -e "\nWait ${waitingTime} hours before re-submit (elapsed ${timeDiff})" >> $outFilename 2>&1
	    exit
	fi
    fi
fi

date > $outFilename 2>&1

###############################
# Mail user in case of problems
###############################
function alertUser()
{
    if [ "$userMail" = "" ]; then
	return
    fi
    alertCode=${2:-"alert!"}
    echo -e "\n$1" >> $outFilename 2>&1
    mail -s "$scriptName: $alertCode" ${userMail} < $outFilename 
}

#########################
# Checks before launching
#########################

# Check that output directory exists
if [ ! -e "${outDir}" ]; then
    alertUser "Output dir ${outDir} not found. Nothing done!"
    exit
fi

# Check alien token
isValidToken=`alien-token-info | grep -c "Token is still valid"`
if [ $isValidToken -eq 0 ]; then
    alertUser "No valid token found. Nothing done!"
    exit
fi
# Check grid proxy
proxyValidity=`xrdgsiproxy info 2>&1 | grep "time left" | cut -d " " -f 6`
if [[ $proxyValidity == "" || $proxyValidity == "0h:0m:0s" ]]; then
    echo  >> $outFilename 2>&1
    alertUser "No valid proxy found. Nothing done!"
    exit
fi

echo -e "\nValid token $isValidToken - proxy $proxyValidity" >> $outFilename 2>&1

cd $outDir
echo -e "\nCurrent directory" >> $outFilename 2>&1
pwd >> $outFilename 2>&1

##########################################
# Get run list from Monalisa (if required)
##########################################
if [ "${prodUrlInMonalisa}" ]; then
#    recoList=`curl -s ${prodUrlInMonalisa} | grep "run numbers" | cut -d "'" -f 4`
#    recoList=${recoList//","/""}

    # Get run list reconstructed with the selected "alirootInMonalisa"
recoList=`curl -s ${prodUrlInMonalisa} | awk -v aliRev="${alirootInMonalisa}" -v qaPassLisa="${qaPassInMonalisa}" ' {
       if ( index($0,"runview")) {
         currLine=$0; gsub("\"","",currLine); split(currLine,arr,"=");
         currRun=arr[3];
         while ( index($0,"table_row") == 0 ) getline;
         for ( iline=0; iline<3; iline++ ) getline;
	 split($0,arr,">"); split(arr[2],finArr,"<");
	 filteredEvents=finArr[1];
         getline;
         while ( index($0,"table_row") == 0 ) getline;
	 for ( iline=0; iline<1; iline++ ) getline;
         alirootRevLine=$0;
	 getline;
	 outDirLisa=$0;
         #if (index(alirootRevLine,aliRev) && filteredEvents>0){
	 if (index(alirootRevLine,aliRev) && index(outDirLisa,qaPassLisa)) print currRun; } } ' | xargs`
#     recoList=`curl -s ${prodUrlInMonalisa} | awk -v aliRev="${alirootInMonalisa}" ' {
#        if ( index($0,"openLive")) {
#          for ( iline=0; iline<4; iline++ ) {
#            getline;
#          }
#          currRun=$0; gsub("\t","",currRun);
#          for ( iline=0; iline<6; iline++ ) {
#            getline;
#          }
# 	 split($0,arr,">"); split(arr[2],finArr,"<");
# 	 filteredEvents=finArr[1];
# 	 for ( iline=0; iline<2; iline++ ) {
#              getline;
#          }
#          alirootRevLine=$0;
#          #if (index(alirootRevLine,aliRev) && filteredEvents>0){
# 	 if (index(alirootRevLine,aliRev)){
#            print currRun; }} } ' | xargs`
    if [ "${recoList}" ]; then
        # run list is not empty
	echo -e "\nBuilding run list form ${prodUrlInMonalisa}" >> $outFilename 2>&1

	# Create temporary run list
	tmpRunListName="tmp_${runListName}"
	for irun in $recoList; do
	    # Remove duplicated runs (don't know why they are present!)
	    if [ -e "${tmpRunListName}" ]; then
		isRunAlreadyThere=`grep -c $irun ${tmpRunListName}`
		if [ $isRunAlreadyThere -gt 0 ]; then
		    continue
		fi
	    fi
	    if [ -e "${excludeRunListName}" ]; then
		isRunExcluded=`grep -c ${irun} ${excludeRunListName}`
		if [ $isRunExcluded -gt 0 ]; then
		    continue
		fi
	    fi
	    echo ${irun} >> ${tmpRunListName}
	done

	# Check if a run list was already there
	if [ -e "${runListName}" ]; then
	    isDiff=`diff -q ${tmpRunListName} ${runListName}`
	    if [ ! "${isDiff}" ]; then
		# If current run list is the same as the previous, do nothing
		echo -e "\nRun list not modified" >> $outFilename 2>&1	
		rm ${tmpRunListName}
		exit
	    fi
	    addRunList=`diff ${tmpRunListName} ${runListName} | grep "<" | xargs`
	    removedRunList=`diff ${tmpRunListName} ${runListName} | grep ">" | xargs`
	    addRunList=${addRunList//"<"/""}
	    removedRunList=${removedRunList//">"/""}
	    echo -e "\nAdded runs: ${addRunList}" >> $outFilename 2>&1
	    echo "Removed runs: ${removedRunList}" >> $outFilename 2>&1
	fi
	mv ${tmpRunListName} ${runListName}
    fi
fi


#######################################
# Find blocking runs and eliminate them.
# We assume that if the last line of the log file is 
# "Accessing file alien://..."
# it means that the file has problem:
# remove the correspondig run and re-launch
#######################################
if [ -e "logMerge.txt" ]; then
    blockingRunStr=`tail -n1 logMerge.txt -n 1`
    # Check that last line contains "alien://"
    isAlienFile=`echo "${blockingRunStr}" | grep -c "alien"`
    if [ $isAlienFile -gt 0 ]; then
	blockingRunStr=`echo ${blockingRunStr} | cut -d "/" -f8`
	blockingRun=${blockingRunStr//000/""}
        # Check that it contains a run number
	if [ "${blockingRun}" != "${blockingRunStr}" ]; then
	    # Kill previous running runCronQA.sh
	    parentProcess=`ps x -o "%r %c " | grep $scriptName | xargs | cut -d " " -f 1`
	    if [ "$parentProcess" ]; then
		childList=`ps x -o "%p %r" | grep "\ ${parentProcess}" | cut -d " " -f1 | xargs`
		kill $childList
		echo "Removing blocking run ${blockingRun}" >> $outFilename 2>&1
		origRunListName="orig_${runListName}"
		if [ ! -e ${origRunListName} ]; then
		    cp -pi ${runListName} ${origRunListName}
		fi
		tmpRunList="tmp_${runListName}"
		cat ${runListName} | grep -v ${blockingRun} > ${tmpRunList}
		mv ${tmpRunList} ${runListName}
		diff ${origRunListName} ${runListName} > "removed_${runListName}"
		alertUser "Removed run ${blockingRun}"
	    fi
	fi
    fi
fi

########
# Run QA
########
runQAopt="-f"
if [ -e "${triggerListName}" ]; then
    runQAopt="${runQAopt} -i ${triggerListName}"
fi

echo -e "\nRunning command:" >> $outFilename 2>&1
echo "  $ALICE_ROOT/PWGPP/MUON/lite/runQA.sh -o ${qaFileName} $runQAopt $runListName ${qaPass} ${alienBaseDir}" >> $outFilename 2>&1

#$ALICE_ROOT/PWGPP/MUON/lite/runQA.sh -o ${qaFileName} $runQAopt $runListName "${qaPass}" "${alienBaseDir}" >> $outFilename 2>&1
$ALICE_ROOT/PWGPP/MUON/lite/runQA.sh -m -a -t -o ${qaFileName} $runQAopt $runListName "${qaPass}" "${alienBaseDir}" >> $outFilename 2>&1


#################
# Check log files
#################
logOut=`find . -name "log*.txt"`
for ifile in $logOut; do
    foundProblem=`grep -c "There was a crash." $ifile`
    if [ $foundProblem -gt 0 ]; then
	alertUser "Found crash in $ifile"
	exit
    fi
done

######################
# Copy output on alien
######################
source /tmp/gclient_env_$UID
#pdfOut=`find . -name "*.pdf" | xargs`
pdfOut=`ls *.pdf terminateRuns/*.pdf | xargs`
for ifile in $pdfOut; do
    currFile=`basename ${ifile}`
    echo -e "\nalien_cp -n $ifile ${outAlienDir}/${currFile}@${alienCloseSE}" >> $outFilename 2>&1
    alien_rm "${outAlienDir/alien:\/\//}/${currFile}"
    alien_cp -n "$ifile" "${outAlienDir}/${currFile}@${alienCloseSE}"
    alien_rm "${outAlienDir/alien:\/\//}/${currFile//.pdf/.root}"
    alien_cp -n "${ifile//.pdf/.root}" "${outAlienDir}/${currFile//.pdf/.root}@${alienCloseSE}"
done

echo -e "\nCheckTime $checkTime" >> $outFilename 2>&1

alertUser "Output written in ${outAlienDir}" "notification"

echo "" >> $outFilename 2>&1
date >> $outFilename 2>&1



# Install crontab: example
# crontab -e
## min(0-59) hour(0-23) day(1-31) month(1-12) dayOfWeek(0-6) command
#15 * 6-7 10 * /users/aliced/stocco/macros/gridAnalysis/muonQA/runCronQA.sh "QAresults_outer.root" > /dev/null 2>&1
#15 * 6-7 10 * /users/aliced/stocco/macros/gridAnalysis/muonQA/runCronQA.sh "QAresults_barrel.root" > /dev/null 2>&1