#!/bin/sh

#usage:
#$1 = base directory where to look for the logs (e.g. /alice/data/2010/LHC10h/000137161/ESDs/pass1)

if [ $# = 0 ]; then
  echo "give the base directory where to look for the logs in arguments"
  exit 3
fi

source /tmp/gclient_env_$UID

if [[ -d "CheckLogs" ]]; then
  echo "a directory containing checks (CheckLogs/) already exist. Please move or remove it"
  exit 3
fi

mkdir CheckLogs
cd CheckLogs

num=0
nEv=0

for dir in `alien_ls $1/`; do
  
  let num+=1
  
# if [ $num -lt 2990 ]
# then
#   continue
# fi
  
  echo "checking file number $num (in $1/$dir)"
  alien_cp alien://$1/$dir/rec.log tmp.txt
  let "nEv += `cat tmp.txt | grep 'type LowMultiplicity ===' | wc -l`"
  echo $nEv > nEvents.txt
  cat tmp.txt | grep "E-AliMUON" >> error.txt
  cat tmp.txt | grep "W-AliMUON" >> warning.txt
  cat tmp.txt | grep "shower event" >> triggerPb.txt
  cat tmp.txt | grep "Too many track candidates" >> trackerPb.txt
  cat tmp.txt | grep "Too many local maxima" >> clusterPb.txt
#  gbbox cat $1/$dir/stdout > tmp.txt 2>&1
#  cat tmp.txt | grep "E-TFileMerger" >> error.txt
  rm -f tmp.txt
  
done

cd ..
