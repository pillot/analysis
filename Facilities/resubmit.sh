#!/bin/bash

# add this line to the crontab (opened with crontab -e) to run it automatically:
#
# 30 * 1-2 4 * ~/Work/Alice/Work/Macros/resubmit.sh #firstJobToConsider > ~/Work/ResubitAutoLogs/`date | sed s/\ /_/g | sed s/:/./g` 2>&1
#
# description of this command:
#
# 30    --> the 30th minutes of the hour
# *     --> every hours
# 1-2 4 --> between 1st and 2nd of April
# *     --> every day within the above limit

if [ $# = 0 ]; then
  echo "you must give the first run number to consider as argument"
  exit
fi

source /Users/pillot/Work/bin/alienenv.sh /Users/pillot/Work/Alice/Work

isValidToken=`alien-token-info | grep -c "Token is still valid"`
if [ $isValidToken -eq 0 ]; then
    echo "No valid token found. Nothing done!"
    exit
fi

proxyValidity=`xrdgsiproxy info 2>&1 | grep "time left" | cut -d " " -f 6`
if [[ $proxyValidity == "" || $proxyValidity == "00h:00m:00s" ]]; then
    echo "No valid proxy found. Nothing done!"
    exit
fi

source $ALICE/Macros/Facilities/gridCommands.sh

gridFindAllFailed $1 resubmit

