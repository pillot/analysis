#!/bin/bash

source /tmp/gclient_env_$UID

###########################
# Get list of master jobs #
###########################
function getMasterList(){
    masterList=`alien_ps -M -b | cut -d " " -f 5 | xargs`;
    eval "$1=\$masterList"
}


##############################################
# Try to determine the first job of a series #
##############################################
function guessFirstJob(){
    echo "Guessing first job"
    jobList="$1"
    invertList=""
    for jobId in $jobList; do
	if [ "$invertList" ]; then
	    invertList="$jobId $invertList"
	else
	    invertList="$jobId"
	fi
    done

    previousJob=0;
    guessedJob=0;
    for jobId in $invertList; do
#	jobIdNum1=${jobId:3}
#	jobIdNum=${jobIdNum1##"0"}
	jobIdNum=${jobId##"0"}
	let "currId = $jobIdNum + 3000"
	if [ $previousJob -gt $currId ]; then
	    break;
	fi
	previousJob=$jobIdNum
	guessedJob=$jobId
    done

    eval "$2=\$guessedJob"
}

###############################################
# Resubmit subjobs in the specified masterjob #
###############################################
function resubmitJob(){
    if [ ! "$1" -o ! "$2" ]; then
	echo "Usage: resubmitJob <masterJob> <errorType>"
	return;
    fi    
    gbbox masterJob $1 -status $2 resubmit 
}


###########################################
# Kill subjobs in the specified masterjob #
###########################################
function killJob(){
    if [ ! "$1" ]; then
	echo "Usage: resubmitJob <masterJob> <errorType>"
	return;
    fi    
    if [ ! "$2" ]; then
      gbbox masterJob $1 kill
      gbbox kill $1
    else  
      gbbox masterJob $1 -status $2 kill
    fi
}


##################################################
# Find masterjob of jobs in the specified status #
##################################################
# function gridFindMaster(){
#     jobStatus="ERROR_V";
#     if [ $1 ]; then
# 	jobStatus="$1";
#     fi

#     minJob=$2;
#     masterList=""
#     getMasterList masterList
#     if [ ! $2 ]; then
# 	guessFirstJob $masterList minJob
# 	#echo "usage: gridFindMaster <minJobNum> [jobStatus]"
# 	#return
#     fi

#     echo "Searching for jobs > $minJob with status $jobStatus"
    
#     jobList=`alien_ps -b -f "$jobStatus" | cut -d "-" -f 2 | cut -d " " -f 1 | xargs`;
#     masterListWithFails=""
#     oldMasterNum=0
#     for ijob in $jobList; do
# 	if ! [[ "$ijob" =~ ^[0-9]+$ ]] ; then
# 	    continue;
# 	fi
# 	if [ $ijob -lt $minJob ]; then
# 	    continue;
# 	fi
# 	masterId=`alien_ps -jdl $ijob | grep MasterJobId`;
# 	masterNum=`echo $masterId | cut -d "\"" -f 2`;
# 	echo "$ijob -> master $masterNum";
# 	if [ $masterNum -ne $oldMasterNum ]; then
# 	    oldMasterNum=$masterNum
# 	    masterListWithFails="${masterListWithFails} $masterNum"
# 	fi
#     done
#     echo "Failing masterjobs:"
#     echo "$masterListWithFails"
#     if [ ! "$masterListWithFails" ]; then
# 	return
#     fi
#     echo ""
#     echo "Resubmit failing jobs? [y/n]"
#     read decision;
#     if [ $decision == "y" ]; then
# 	for imaster in $masterListWithFails; do
# 	    resubmitJob $imaster $jobStatus
# 	done
#     fi
# }


###################################################
# Print status of the masterjobs with an id >= $1 #
###################################################
function gridCheckMasters(){
    if [[ ! $1 ]]; then
        echo "give the id of the first masterjob to be consider"
	return
    fi
    masterList=""
    getMasterList masterList
    for ijob in $masterList; do
	if [ $ijob -lt $1 ]; then
	    continue;
	fi
	masterstatus=`alien_ps -b -M -j $ijob | tr -s " " | cut -d" " -f5`
	subjobstatus=`alien_ps -b -m $ijob | tr -s " " | cut -d" " -f5 | xargs`
	nsubjobs=0
	if [[ $subjobstatus ]]; then
	    nsubjobs=`echo -e ${subjobstatus//" "/"\n"} | wc -l | xargs`
	fi
	ndone=`echo -e ${subjobstatus//" "/"\n"} | grep -c ^D`
	nwaiting=`echo -e ${subjobstatus//" "/"\n"} | grep -c ^W`
	nerror=`echo -e ${subjobstatus//" "/"\n"} | grep -c ^E`
	nzombie=`echo -e ${subjobstatus//" "/"\n"} | grep -c ^Z`
	let "nother=$nsubjobs - $ndone - $nwaiting - $nerror - $nzombie"
	echo "masterjob $ijob in status $masterstatus -- subjobs:$nsubjobs (D:$ndone, W:$nwaiting, E:$nerror, Z:$nzombie, O:$nother)"
    done
}


###########################################################
# Find and resubmit all jobs in error $1 with an id >= $2 #
# If $3 == "resubmit" then automatically resubmit them    #
# Else ask before resubmission                            #
###########################################################
function gridFindFailed(){
    jobStatus="ERROR_V";
    if [ $1 ]; then
	jobStatus=$(echo $1 | tr [:lower:] [:upper:]);
    fi

    minJob=$2;

    masterList=""
    if [[ ! $2 ]]; then
	getMasterList masterList
	guessFirstJob "$masterList" minJob
    fi

    declare -a errStatusFull=("ERROR_V" "ERROR_E" "ERROR_IB" "ERROR_SV" "EXPIRED" "ZOMBIE" "ERROR_VN" "ERROR_RE" "ERROR_VT" "ERROR_A" "ERROR_SPLT")
    declare -a errStatusShort=("EV" "EE" "EIB" "ESV" "EXPIRED" "Z" "EVN" "ERE" "EVT" "EA" "ESPLT")

    nStatus=${#errStatusFull[@]}
    jobStatusShort=""
    for (( ist=0; ist<$nStatus; ist++ )){
	if [ ${errStatusFull[$ist]} == "$jobStatus" ]; then
	    jobStatusShort="${errStatusShort[$ist]}"
	    break
	fi
    }

    echo "Searching for jobs > $minJob with status $jobStatus ($jobStatusShort)"

    masterListWithFails="";
    masterListFails="";
    # method 1: more efficient when few jobs failing
    jobList=`alien_ps -b -l 10000 -f "$jobStatus" | sed s/"  *"/" "/g | cut -d " " -f 3 | xargs`;
    #jobList=`gbbox 'ps -AS' | grep " $jobStatusShort" | sed s/"  *"/" "/g | cut -d " " -f 2 | xargs`;
    oldMasterNum=0
    for ijob in $jobList; do
	if [[ "$ijob" =~ ^[0-9]+$ ]] ; then
	    if [ $ijob -lt $minJob ]; then
		continue;
	    fi
	    numSubJob=`alien_ps -b -m $ijob | wc -l`
	    if [[ $numSubJob -eq 0 ]] ; then
#	        if [[ $masterNum -gt 169958885 ]]; then
#		    continue;
#	        fi
		masterListFails="${masterListFails} $ijob"
	    fi
	else
	    ijob=${ijob:1}
	    if [ $ijob -lt $minJob ]; then
		continue;
	    fi
	    masterId=`alien_ps -jdl $ijob | grep MasterJobId`;
	    masterNum=`echo $masterId | cut -d "\"" -f 2`;
	    if [[ ! -n "$masterNum" ]] || [[ $masterNum -lt $minJob ]]; then
		continue;
	    fi
#	    if [[ $masterNum -gt 169958885 ]]; then
#		continue;
#	    fi
	    echo "$ijob -> master $masterNum";
	    oldMasterNum=`echo $masterListWithFails | grep $masterNum`
	    if ! [[ -n "$oldMasterNum" ]]; then
		masterListWithFails="${masterListWithFails} $masterNum"
	    fi
	fi
    done
    if [ "$masterListWithFails" ]; then
        echo "--> Failing masterjobs with subjobs: $masterListWithFails"
	decision="y"
	if ! [[ $3 == "resubmit" ]]; then
          echo "Resubmit failing jobs? [y/n]"
          read decision;
	fi
        if [ $decision == "y" ]; then
	    for imaster in $masterListWithFails; do
	        resubmitJob $imaster $jobStatus
#                killJob $imaster $jobStatus
#                killJob $imaster
	    done
        fi
    fi
    if [ "$masterListFails" ]; then
        echo "--> Failing masterjobs without subjobs: $masterListFails"
	decision="y"
	if ! [[ $3 == "resubmit" ]]; then
          echo "Resubmit failing jobs? [y/n]"
          read decision;
	fi
        if [ $decision == "y" ]; then
            for imaster in $masterListFails; do
                gbbox resubmit $imaster
#                killJob $imaster
	    done
        fi
    fi
}


###########################################################
# Find and resubmit all jobs in error $1 with an id >= $2 #
# If $3 == "resubmit" then automatically resubmit them    #
# Else ask before resubmission                            #
###########################################################
function gridFindFailed2(){
    jobStatus="ERROR_V";
    if [ $1 ]; then
	jobStatus=$(echo $1 | tr [:lower:] [:upper:]);
    fi

    minJob=$2;

    masterList=""
    if [[ ! $2 ]]; then
	getMasterList masterList
	guessFirstJob "$masterList" minJob
    fi

    declare -a errStatusFull=("ERROR" "ERROR_V" "ERROR_E" "ERROR_IB" "ERROR_SV" "EXPIRED" "ZOMBIE" "ERROR_VN" "ERROR_RE" "ERROR_VT" "ERROR_A" "ERROR_SPLT")
    declare -a errStatusShort=("E" "EV" "EE" "EIB" "ESV" "EXPIRED" "Z" "EVN" "ERE" "EVT" "EA" "ESPLT")

    nStatus=${#errStatusFull[@]}
    jobStatusShort=""
    for (( ist=0; ist<$nStatus; ist++ )){
	if [ ${errStatusFull[$ist]} == "$jobStatus" ]; then
	    jobStatusShort="${errStatusShort[$ist]}"
	    break
	fi
    }

    echo "Searching for jobs > $minJob with status $jobStatus ($jobStatusShort)"

    listFails="";
    # method 1: more efficient when few jobs failing
    jobList=`gbbox 'ps -AS' | grep " $jobStatusShort" | sed s/"  *"/" "/g | cut -d " " -f 2 | xargs`;
    oldMasterNum=0
    for ijob in $jobList; do
	if [[ "$ijob" =~ ^[0-9]+$ ]] ; then
	    if [ $ijob -lt $minJob ]; then
		continue;
	    fi
	    numSubJob=`gbbox masterJob $ijob | grep "In total, there are" | cut -d " " -f 5`
	    if [[ $numSubJob -eq 0 ]] ; then
#	        if [[ $masterNum -gt 169958885 ]]; then
#		    continue;
#	        fi
		listFails="${listFails} $ijob"
	    fi
	else
	    ijob=${ijob:1}
	    if [ $ijob -lt $minJob ]; then
		continue;
	    fi
	    masterId=`alien_ps -trace $ijob | grep "Master Job is" | cut -d "[" -f 4 | cut -d " " -f 4 | sed s/"\]"//g`
	    if [ $masterId -lt $minJob ]; then
		continue;
	    fi
	    listFails="${listFails} $ijob"
	fi
    done
    if [ "$listFails" ]; then
        echo "--> Failing jobs: $listFails"
	decision="y"
	if ! [[ $3 == "resubmit" ]]; then
          echo "Resubmit failing jobs? [y/n]"
          read decision;
	fi
        if [ $decision == "y" ]; then
            for iJob in $listFails; do
                gbbox resubmit $iJob
#                killJob $imaster
	    done
        fi
    fi
}


########################################################
# Find and resubmit all jobs in error with an id >= $1 #
# If $2 == "resubmit" then automatically resubmit them #
# Else ask before resubmission                         #
########################################################
function gridFindAllFailed(){
    declare -a errStatusFull=("ERROR_V" "ERROR_E" "ERROR_IB" "ERROR_SV" "EXPIRED" "ZOMBIE" "ERROR_VN" "ERROR_RE" "ERROR_VT" "ERROR_A" "ERROR_SPLT")
    
    nStatus=${#errStatusFull[@]}
    for (( ist=0; ist<$nStatus; ist++ )){
	echo "------------------------------------"
        gridFindFailed ${errStatusFull[$ist]} $1 $2
        echo ""
    }
}


########################################################
# Find and resubmit all jobs in error with an id >= $1 #
# If $2 == "resubmit" then automatically resubmit them #
# Else ask before resubmission                         #
########################################################
function gridFindAllFailed2(){
    echo "------------------------------------"
    gridFindFailed2 "ERROR" $1 $2
    echo ""
}


####################
# Find worker node #
####################
function gridFindWN(){
    if [ ! $1 ]; then
	echo "Usage: findWN <masterjob number> [psOptions]"
	return
    fi
    psOption="-f ERROR_%"
    if [ "$2" ]; then
	psOption="$2"
    fi
    alien_ps -m $1 $psOption -F l
#     subRunList=`alien_ps -b -m $1 $psOption | cut -d " " -f 5`
#     for jobId in $subRunList; do 
# 	jobNum=${jobId:1}
# 	echo $jobNum
# 	psOut=`alien_ps -trace $jobNum all | grep "transition to STARTED"`
# 	echo $psOut;
#     done
}


####################
# Count chunks     #
####################
function gridCountChunks(){
    if [ ! $1 ]; then
	echo "Usage: gridCountChunks <masterjob number> [psOptions]"
	return
    fi
    psOption="-f DONE"
    if [ "$2" ]; then
	psOption="$2"
    fi
    showDetails=0
    if [ $3 ]; then
	showDetails=1;
    fi
    nKilledJobs=`alien_ps -trace $1 | grep -c "killing"`
    nChunksInXml="?"
    if [ $4 ]; then
	collectionFile=`alien_ps -trace $1 all | grep xml | cut -d ":" -f 5 | cut -d "," -f 1`
	localTmpXml="/tmp/tmpCollectionName.txt"
	alien_cp -s "alien://$collectionFile" "/tmp/tmpCollectionName.txt"
	nChunksInXml=`grep -c "event name" $localTmpXml`;
	rm $localTmpXml
    fi
    subRunList=`alien_ps -b -m $1 $psOption | cut -d " " -f 5`
    totalPerRun=0;
    for jobId in $subRunList; do 
	jobNum=${jobId:1}
	nChunks=`alien_ps -jdl $jobNum | grep -c "LF:/alice/data"`
# 	lfList=`alien_ps -jdl $jobNum | grep "LF:/alice/data"`
#  	nChunks=`echo $lfList | grep -c "ESD"`
# 	if [ $nChunks -eq 0 ]; then
# 	    nChunks=`echo $lfList | grep -c "AOD"`  
# 	fi
	if [ $showDetails -eq 0 ]; then
	    echo "job $jobNum   chunks $nChunks"
	fi
	let "totalPerRun += $nChunks"
    done
    echo "Master $1: chunks analyzed $totalPerRun  total $nChunksInXml  killed jobs $nKilledJobs";
}


####################
# Chunk summary    #
####################
function gridChunkSummary(){
    psOption=$1
    minJob=$2;
    if [ ! $2 ]; then
	guessFirstJob minJob
    fi
    echo "Searching for jobs > $minJob with options $psOption"
    jobList=`alien_ps -M -b | cut -d " " -f 5 | xargs`;
    for jobId in $jobList; do
	if [ $jobId -lt $minJob ]; then
	    continue;
	fi
	gridCountChunks "$jobId" "$psOption" 0
    done
}


# ###########################
# # Register uncopied SAVED #
# ###########################
# function gridRegisterSaved(){
#     #jobList=`alien_ps -b -f "SAVED" | cut -d "-" -f 2 | cut -d " " -f 1 | xargs`;
#     jobList="47516210" # REMEMBER TO CHANGE
#     for jobId in $jobList; do
# 	outDir=`alien_ps -jdl $jobId | grep "OutputDir =" | cut -d "\"" -f 2`
# 	alien_mkdir -p $outDir
# 	alien_cd $outDir
# 	alien_registerOutput $jobId
#     done
# }

##############################
# Copy files to specified SE #
##############################
function copyToSE() {
    if [[ ! "$1" || ! "$2" ]]; then
	echo "Usage: copyToSE <path> <pattern> [destinationSE]"
    fi
    filePath="$1"
    filePattern="$2"

    finalSE="$alien_CLOSE_SE"
    if [ $3 ]; then
	finalSE="$3"
    fi
    if [ ! "$finalSE" ]; then
	echo "Cannot find SE. Please define it explicitely"
	return
    fi
    echo "Final storage element $finalSE"    
    
    fileList=`alien_find "${filePath}" "${filePattern}" | grep -v "files found"`

    for file in $fileList; do
	#isInSE=`alien_whereis $file 2>&1 | grep -c "$finalSE"`
	#if [ $isInSE -gt 0 ]; then
	    continue
	#fi
	#alien_cp "alien://${file}" "alien://${file}@${finalSE}"
	    echo alien_mirror "${file}" "${finalSE}"
    done
}


###############################
# Remove chunks after merging #
###############################
function gridCleanMerged() {
    if [ ! $1 ]; then
	echo "Usage: gridCleanMerged <alien-outdir>"
	return
    fi

    baseOutDir=$1;

    mergedFileKey="*Stage*.zip"

    mergedList=`alien_find ${baseOutDir} ${mergedFileKey} | grep -v "files found"`
    for file in $mergedList; do
	currDir=`dirname ${file}`
	#subDirList=`alien_find ${currDir}/ 0* | grep -v "files found"`
	#for subDir in $subDirList; do
	#echo "alien_rm_finto ${subDir}"
	#done
	rmDir="${currDir}/0*"
	echo "Remove directories ${rmDir}? [y/n]"
	read decision
	if [ "${decision}" == "y" ]; then
	    alien_rmdir "${rmDir}"
	fi
    done
}

###############################
# Remove chunks after merging #
###############################
function gridMoveMerged() {
    if [ ! $1 ]; then
	echo "Usage: gridMoveMerged <alien-outdir>"
	return
    fi
    baseOutDir=$1;

    outFileName="AnalysisResults.root"
    if [ $2 ]; then
	outFileName="$2"
    fi

    outRunList=`alien_ls -F ${baseOutDir} | grep "/"`
    outRunList=${outRunList//\//""}

    for irun in $outRunList; do
	mergedFile="${baseOutDir}/${irun}/${outFileName}"
	alien_mv ${mergedFile} ${mergedFile//".root"/"-merged.root"}
    done
}


###############################
# Remove chunks after merging #
###############################
function gridCheckMerged() {
    if [ ! $1 ]; then
	echo "Usage: gridCleanMerged <alien-outdir>"
	return
    fi
    baseOutDir=$1;

    outFileName="AnalysisResults.root"
    if [ $2 ]; then
	outFileName="$2"
    fi

    outFileKey="root_archive.zip"
    mergedFileKey="*Stage*.zip"


    #outRunList=`alien_ls $baseOutDir | grep 000`
    outRunList=`alien_ls -F ${baseOutDir} | grep "/"`
    outRunList=${outRunList//\//""}
    mergedList=`alien_find ${baseOutDir} ${mergedFileKey} | grep -v "files found"`
    allOutList=`alien_find ${baseOutDir} ${outFileKey} | grep -v "files found"`

    runToMerge=""
    incompleteRuns=""

    currList="$allOutList"
    currListAux=""

    yesToAll=0

    # Loop on runs in output
    for irun in $outRunList; do
	echo "Checking run ${irun}"
	nMerged=`echo ${mergedList} | grep -c $irun`
	# Check if run was merged
	if [ $nMerged -gt 0 ]; then

	    # Check if all sub-jobs were produced
	    # (Warning: it checks for missing chuncks,
	    # but if the last are missing it cannot find it
	    nOuts=0
	    subList=""
	    for file in $currList; do
		isMatch=`echo $file | grep -c $irun`
		if [ $isMatch -gt 0 ]; then
		    let "nOuts++";
		    subList="${subList} ${file}"
		    else
		    currListAux="${currListAux} ${file}"
		fi
	    done
	    currList="${currListAux}"
	    currListAux=""
	    for ((isub=1; isub<=$nOuts; isub++)); do
		currSub="${isub}"
		if [ $isub -lt 10 ]; then
		    currSub="00${isub}"
		elif [ $isub -lt 100 ]; then
		    currSub="0${isub}"
		fi
		nMatch=`echo ${subList} | grep -c "${irun}/${currSub}"`
		if [ $nMatch -eq 0 ]; then
		    incompleteRuns="${incompleteRuns} ${irun}"
		    break
		fi
	    done

	    mergedFile=`alien_find ${baseOutDir} ${irun}/${outFileName} | grep -v "files found"`
	    if [ "${mergedFile}" ]; then
		if [ ${yesToAll} -eq 0 ]; then
		    echo "Move ${mergedFile}? [y/n/a]"
		    read decision
		fi
		if [ "${decision}" == "a" ]; then
		    yesToAll=1;
		    decision="y"
		fi
		if [ "${decision}" == "y" ]; then
		    alien_mv ${mergedFile} ${mergedFile//".root"/"-merged.root"}
		fi
	    fi

	    # If all chuncks were produced and merged
	    # Remove the chunks
# 	    if [ $nMatch -gt 0 ]; then
# 		rmDir="${baseOutDir}/${irun}/0*"
# 		echo "Remove directories ${rmDir}? [y/n]"
# 		read decision
# 		if [ "${decision}" == "y" ]; then
# 		    echo alien_rmdir "${rmDir}"
# 		fi
# 	    fi
	else
	    runToMerge="${runToMerge} $irun"
	fi
    done

    echo "Runs to merge: $runToMerge"
    echo "Incomplete runs: $incompleteRuns"

}