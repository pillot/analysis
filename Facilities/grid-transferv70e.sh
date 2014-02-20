#!/bin/bash

# Release notes v53
# copyerrorfiles: copy from next SE is not started before all existing copy
# processes are finished.
# watchdog: sleeptime as parameters in seconds in config file
# copyerrorfiles: wrongly copied files are removed from Lustre
# automatic proxy renewal included
# independet of number of replicas at least copy 2 times
# closese not copied to .alienshrc anymore (kept in memory instead"
# remove wrongly copied files activated
# v60: SE list based on more than 2 files from $errorfiles
# v61: blacklisted SE can be specified in config file, will be excluded from SElist
# v62: unzip check right after file has been copied
# v63: improved se list creation
# v64: bug removed in errorcopy procedure (main)
# v65: but removed in createselist, se based on 3 files only
#    : fixed also bug in main: createerrorfile setting of aliencloseSE
# v65a-nozip: no test for zip files, enhanced directory flexibility in copyfile and copyerrorfile
# v65b-nozip: config file and errorfiles are copied to the processed dir after everything is done
# v65c-nozip: introduced else if in line 150 to fix unary operator problem
# v65d-nozip: introduced check if all backup error files are existing
# v65e-nozip: skipping already copied files (see v68)
# v65f-nozip: config file and errorfiles are moved to the processed dir after everything is done
# v65g-nozip: date output included in log file (watchdog function)
# v65h-nozip: chgrp and chmod added after mkdir to ensure proper permissions
# v65j-nozip: umask introduced. chgrp after each mkdir and after file check
# v65k-nozip: catch empty lines in config file 
# v65l-nozip: new copy mode (12 random files from a run) - config file extended
# v65m-nozip: replacing alien -exec find by alien_find
# v65n-nozip: changed to new transfer infrastructure and AliEn v218
# v65o-nozip: changed back alien_find by alien -exec find due to errors with alien_find
# v65p-nozip: introduced x3 in copyfile function to enable longer directory names
# v65pq-nozip: no proxy-renewal anymore after each run. Is done outside of the script in a cronjob
# v66r: introduced check if file pattern is zip file. Unzip check reactivated
#       introduced function waitfortermination
#       introduced usage of waitfortermination function in code
#       copyerrorfiles: improved checking of existence of errorfiles
#       reactivated that errorfilesbackup is copied back to errorfiles ofter processing
#       activated automatic removal of wrongly copied files
# v67r: enabled list mode copy as option via config file
# v67s: proxyrenewal also in listmode
# v68s: proxyrenewal also in additional listmode copy. 
#       proxyrenewal only sources the Token Environment
#       removeerrorfiles activated, but in function actual removal is not yet done.
# v69: if no runs are in the config file or if the token creation does
#      not work error messages are produced. No switch to list mode is done. The program exists.
#      exit code 1: no or wrong transfer mode
#      exit code 2: no token
# v69b: but removed if [[ "$copymode" != "$listmode" ]
# v70: apply summary check to see how many files have been copied
#      apply "" to comparison statements in lines 231, 236, 240
#      listmode copy only when copy mode listmode and errorfile exists  
# v70b: random mode via flexible input parameter from config file in percentage from run
# v70c: proxy renewal in listmode for large data sets
# v70d: proxy password as input not needed anymore, removed some commented lines
# v70e: changed to AliEn v219. replace alien -exec find by alien_find

umask u=rwx,g=rwx,o=rx

config=$1

# removing potential blank lines from config file
if [ "$config" != "/lustre/alice/alien/$INPUTDIR/" ] ; then
 mv $config $config"backup"
 grep -v ^$ $config"backup" > $config
 rm $config"backup"
fi 

declare -a array
declare element_count
declare -a current_streams
declare stream_count
declare -a dirstring
declare dirstring_count
declare -a basestring
declare basestring_count
declare -a searray
declare secount
declare -a searray_tmp
declare secount_tmp
declare -a filearray
declare filearray_count

noerror="No errors detected in compressed data of"

# source /misc/kschwarz/.alienshrc

function configuration {
counter=1
while read line || [ "$line" != "" ]
 do
 line[counter]=$line
 counter=$(($counter+1))
done < $config

basedir=${line[1]}
echo "basedir = $basedir"
errorfiles=${line[2]}
echo "errorfiles = $errorfiles"
maxstreams=${line[3]}
echo "maxstreams=$maxstreams"
localbase=${line[4]}
echo "local basedir = $localbase"
filepattern=${line[5]}
echo "file pattern = $filepattern"
# check for zip file
if [ "${filepattern:(-4)}" == ".zip" ] ; then
 echo "it is a zip file. Unzip check activated"
 ziptest=1
 else
 echo "it is not a zip file. Unzip check disabled"
 ziptest=0
fi
sleeptime=${line[6]}
echo "sleep time = $sleeptime"
blackse=${line[7]}
echo "blacklisted se = $blackse"
copymode=${line[8]}
echo "copy mode = $copymode"
randompercent=${line[9]}
echo "random percentage = $randompercent"

for ((a=10; a<$counter; a++))
  do
  run[a-9]=${line[a]}
  echo "configuration: run[$(($a-9))] = ${line[a]}"
done
}

function createfilelist {
echo "createfilelist:run[$a] = ${run[a]}"
echo "createfilelist:searching for $filepattern"
echo "creatfilelist:doing alien_find $basedir/${run[a]} $filepattern"
array=($(alien_find $basedir/${run[a]} $filepattern)) 
element_count=${#array[*]}
echo "createfilelist:element_count = $element_count"
if [ $copymode == "random" ]; then
 echo "createfilelist: random copy mode"
 let randomnumber=$element_count*$randompercent/100
 echo "createfilelist: $randomnumber out of $element_count files will be transferred."
 if [ $element_count -lt $randomnumber ]; then
  randmax=$element_count
 else
  randmax=$randomnumber
 fi
 for ((r=0; r<$randmax; r++))
   do
    number=$RANDOM
    let "number %= $element_count"
    array[r]=${array[number]}
    echo "file number $r = ${array[r]}"
   done
   element_count=$randmax
 else
  echo "createfilelist: standard copy mode"
fi
}

function testmaindir {
if [ -d "$localbase/$basedir/${run[a]}" ]; then
 echo "testmaindir:directory existing"
else
echo "testmaindir:doing mkdir -p $localbase/$basedir/${run[a]}"
mkdir -p "$localbase/$basedir/${run[a]}"
chgrp --quiet alidata "$localbase/$basedir/${run[a]}"
fi
}

function testsubdir {
if [ -d "$localbase/$basedir/${run[a]}/$1" ]; then
 echo "testsubdir:directory existing"
else
echo "testsubdir: doing  mkdir -p $localbase/$basedir/${run[a]}/$1"
mkdir -p "$localbase/$basedir/${run[a]}/$1"
chgrp --quiet alidata "$localbase/$basedir/${run[a]}/$1"
fi
}

function checkfile {
echo "checkfile:checking ${array[$1]}"
cataloguesize=$(alien -exec ls -l ${array[$1]}| awk '{ print $4 }')
localsize=$(ls -l $localbase/${array[$1]}| awk '{ print $5 }')
chgrp --quiet alidata "$localbase/${array[$1]}"
echo "checkfile:cataloguesize = $cataloguesize"
echo "checkfile:local size = $localsize"
echo "checkfile: group changed to alidata"
return=$(unzip -q -t $localbase/${array[$1]})
echo "checkfile: unzip - $return"
if [ "$cataloguesize" != "$localsize" ] ; then
 # do effective unzip check only if file pattern contains zip file
 if [ $ziptest != 1 ]; then
  echo "checkfile: adding ${array[$1]} to $errorfiles"
  echo "${array[$1]}" >> $errorfiles
 else
  if [[ "$cataloguesize" != "$localsize" ]] || [[ $return != $noerror* ]]; then
   echo "checkfile: unzip failed - adding ${array[$1]} to $errorfiles"
   echo "${array[$1]}" >> $errorfiles
  fi 
 fi  
else
 echo "checkfile:HURRAY !!!"
fi        
}

function watchdog {
echo "watchdog: sleeping for $sleeptime seconds"
date
sleep $sleeptime
if [ "$(ps -p $1 | grep $1)" ]; then
 kill -9 $1
 echo "watchdog:killed process $1"
else
 echo "watchdog:process $1 terminated naturally"
fi
checkfile $2
}

function copyfile {
# length=$(echo ${array[c]} | awk '{ print length($1)}')
# x=$(echo ${array[c]} | awk '{ print substr($1,'$(($length-19))',3)}')
dirstring=($(echo ${array[c]} | sed 's/\// /g'))
dirstring_count=${#dirstring[*]}
echo "copyfile:dirstring_count=$dirstring_count"
x=${dirstring[$((dirstring_count-2))]}
x1=${dirstring[$((dirstring_count-3))]}
x2=${dirstring[$((dirstring_count-4))]}
x3=${dirstring[$((dirstring_count-5))]}
echo "copyfile:x=$x, x1=$x1, x2=$x2, x3=$x3"
echo "copyfile:run[$a] = ${run[a]}"
if [ $x == ${run[a]} ]; then
 x=""
 x1=""
 x2=""
 x3=""
 echo "copyfile: x, x1, x2, x3 set to $x, $x1, $x2, $x3"
fi

if [ "$x1" == "${run[a]}" ]; then
 x1=""
 x2=""
 x3=""
 echo "copyfile: x1, x2, x3 set to $x1, $x2, $x3"
elif [ "$x2" == "${run[a]}" ]; then
 x2=""
 x3=""
 echo "copyfile: x2, x3 set to $x2, $x3"
elif [ "$x3" == "${run[a]}" ]; then
 x3=""
 echo "copyfile: x3 set to $x3"
else
 echo "copyfile: no x value changed" 
fi

testsubdir $x3/$x2/$x1/$x
if [ -e "$localbase/${array[c]}" ]; then
 echo "copyfile: file $localbase/${array[c]} exists. Skipping copy process."
else 
 echo "copyfile:copying ${array[c]}"
 echo "copyfile:doing alien_cp -v -t 2 alien:${array[c]} file:/$localbase/$basedir/${run[a]}/$x3/$x2/$x1/$x"
 alien_cp -v -t 2 alien:${array[c]} file:/$localbase/$basedir/${run[a]}/$x3/$x2/$x1/$x &
 watchdog $! $c &
fi 
}

function copyerrorfile {
echo "copyerrorfile: array[$c] = ${array[c]}"
dirstring=($(echo ${array[c]} | sed 's/\// /g'))
dirstring_count=${#dirstring[*]}
echo "copyerrorfile: dirstring_count = $dirstring_count"
filename=${dirstring[$((dirstring_count-1))]}
echo "copyerrorfile: file name = $filename"
slength=$(echo ${#array[c]})
echo "copyerrorfile: string length = $slength"
flength=$(echo ${#filename})
echo "copyerrorfile: length of file name = $flength"
copydir=${array[c]:0:$slength-$flength}
echo "copyerrorfile: directory to be copied to = $copydir"
echo "copyerrorfile: testing subdir"
if [ -d "$copydir" ]; then
 echo "copyerrorfile: sub directory existing"
else
 echo "copyerrorfile: doing  mkdir -p $localbase/$copydir"
 mkdir -p "$localbase/$copydir"
# check first with Jacek 
# chgrp -R --quiet alidata "$localbase/$basedir/$run/$x"
# chmod -R --quiet g+w "$localbase/$basedir/$run/$x"
fi
echo "copyerrorfile:copying ${array[c]}"
echo "copyerrorfile:doing alien_cp -v -m alien:${array[c]} file:/$localbase/$copydir/$filename"
alien_cp -v -m alien:${array[c]} file:/$localbase/$copydir/$filename &
watchdog $! $c &
}

function checkstreams {
current_streams=($(ps -ef | grep alien_cp | awk '{ print $2 }'))
stream_count=${#current_streams[*]}
  echo "checkstreams:stream_count = $stream_count"
   if [ "$stream_count" -gt "$maxstreams" ]; then
    echo "checkstreams:too many streams"
    echo "checkstreams:max streams = $maxstreams"
    echo "checkstreams: sleeping for 1 minute"
    sleep 60
     else
      echo "checkstreams:everything ok"
   fi
}

function waitfortermination {
checkstreams
while [ $stream_count -gt 2 ] ; do
   echo "main: stream_count = $stream_count"
   checkstreams
   sleep 60
done
}

function createselist {
searray_tmp=""
secount=0
secompare1="ALICE::"
secompare2="Alice::"

# reading first line of $errorfiles
# Redirecting stdin using 'exec'.
exec 6<&0          # Link file descriptor #6 with stdin.
                   # Saves stdin.
exec < $errorfiles   # stdin replaced by file "$errorfiles"
read e1            # Reads first line of file "$errorfiles"
read e2            # Reads second line of file "$errorfiles"
read e11            # Reads second line of file "$errorfiles"
#read e21           # Reads second line of file "$errorfiles"
#read e31            # Reads second line of file "$errorfiles"
echo "createselist: SE list based on files $e1, $e2, $e11"
#  echo "createselist: SE list based additionally on files $e21, and $e31"
exec 0<&6 #6<&      # releasing stdin

stderr="$(alien -exec whereis $e1 2>&1)"
searray_tmp=($(echo $stderr))
secount_tmp=${#searray_tmp[*]}
echo "createselist: secount_tmp = $secount_tmp"

for ((a=1; a<$secount_tmp; a++))
 do
#   echo "createselist: searray_tmp[$a] = ${searray_tmp[a]}"
  
   if [[ ${searray_tmp[a]} == $secompare1* ]] || [[ ${searray_tmp[a]} == $secompare2* ]]; then
   echo "createselist: adding ${searray_tmp[a]} to selist"
   secount=$(($secount+1))
   searray[$secount]=${searray_tmp[a]}
    if [ ${searray[$secount]} == $blackse ]; then
      echo "searray[$secount] is blacklisted. Setting to random"
      searray[$secount]="random"
      echo "createselist: SE[$a] now set to ${searray[$secount]}"
    fi  
   fi   
done

searray_tmp=""
stderr="$(alien -exec whereis $e2 2>&1)"
searray_tmp=($(echo $stderr))
secount_tmp=${#searray_tmp[*]}
echo "createselist: secount_tmp (2) = $secount_tmp"

for ((a=1; a<$secount_tmp; a++))
 do
#   echo "createselist: searray_tmp[$a] = ${searray_tmp[a]}"
  
   if [[ ${searray_tmp[a]} == $secompare1* ]] || [[ ${searray_tmp[a]} == $secompare2* ]]; then
   echo "createselist: adding ${searray_tmp[a]} to selist"
   secount=$(($secount+1))
   searray[$secount]=${searray_tmp[a]}
    if [ ${searray[$secount]} == $blackse ]; then
     echo "searray[$secount] is blacklisted. Setting to random"
     searray[$secount]="random"
     echo "createselist: SE[$a] now set to ${searray[$secount]}"
    fi
   fi    
done

searray_tmp=""
stderr="$(alien -exec whereis $e11 2>&1)"
searray_tmp=($(echo $stderr))
secount_tmp=${#searray_tmp[*]}
echo "createselist: secount_tmp (3) = $secount_tmp"

for ((a=1; a<$secount_tmp; a++))
 do
#   echo "createselist: searray_tmp[$a] = ${searray_tmp[a]}"
  
   if [[ ${searray_tmp[a]} == $secompare1* ]] || [[ ${searray_tmp[a]} == $secompare2* ]]; then
   echo "createselist: adding ${searray_tmp[a]} to selist"
   secount=$(($secount+1))
   searray[$secount]=${searray_tmp[a]}
    if [ ${searray[$secount]} == $blackse ]; then
      echo "createselist: searray[$secount] is blacklisted. Setting to random"
      searray[$secount]="random"
      echo "createselist: SE[$a] now set to ${searray[$secount]}"
    fi      
   fi   
done

echo "the following list of SEs has been created"
for ((a=1;a<=$secount;a++))
 do
  echo "SE[$a] = ${searray[$a]}"
 done

}


function removeerrorfiles {
 while read line || [ "$line" != "" ]
  do
   echo "removeerrorfiles: removing $localbase/$line"
   # rm $localbase/$line
  done < $errorfiles
}

function proxyrenewal {
  echo "proxyrenewal: doing proxy renewal"
  echo "AliEnROOT = $ALIEN_ROOT"
  echo "AliEn used is: $(which alien)"
  echo "proxyrenewal: sourcing environment"
  source /tmp/gclient_env_$UID
  tokeninfo=$(alien-token-info)
  if [ "$tokeninfo" == "No Token found!" ]; then
   echo "ERROR: no valid token available !!!"
   exit 2  
  fi   
}

# MAIN PROGRAM

configuration

echo "ALIEN_ROOT = $ALIEN_ROOT"
echo "AliEn Executable = $(which alien)"
echo "gbbox used = $(which gbbox)"

# main copy procedure
# a = counter for run numbers
# if copy mode = listmode then this part is skipped
if [ $copymode != "listmode" ]; then
 for ((a=1; a<$counter-9; a++))
  do
   echo "main: a = $a"
   proxyrenewal
   createfilelist
   testmaindir
   for ((c=0; c<$element_count; c++))
    do
     echo "main: c = $c"
     copyfile
     checkstreams
     checkfivehundred=$(( $c % 500 ))
     if [ $checkfivehundred = 0 ]; then
      echo "c = $c"
      proxyrenewal 
     fi       
   done
 done 
waitfortermination
fi

# wait until old copy processes are terminated
if [ "$copymode" == "" ]; then
 echo "ERROR: no copy mode selected !!!"
 exit 1  
fi 

if [[ "$copymode" != "listmode" ]] && [[ "$copymode" != "standard" ]] && [[ "$copymode" != "random" ]]; then
 echo "ERROR: no valid copy mode selected !!!"
 exit 1
fi

if [ -e "$errorfiles" ] || [ "$copymode" == "listmode" ]; then
 # copy in listmode only when errorfile exists or copymode = listmode
 echo "copy in listmode"
 # try to copy files from error list once more
 # one try per SE
 # source environment
 proxyrenewal

 createselist


 # print selist
 # echo "printing selist"
 errorfilesbackup="$errorfiles-backup0"
 cp $errorfiles $errorfilesbackup
 echo "main: copying $errorfiles to $errorfilesbackup"
 rm $errorfiles
 echo "main: removing old file $errorfiles"


 for ((a=1; a<=$secount; a++))
  do
   echo "main: a=$a, secount=$secount"
   echo "main: copying from SE $a = ${searray[$a]}"
 #  echo "main: writing SE to $HOME/.alienshrc"
 #  echo "export alien_CLOSE_SE=${searray[a]}" > $HOME/.alienshrc
 #  echo "main: sourcing $HOME/.alienshrc"
 #  source $HOME/.alienshrc
   echo "main: export alien_CLOSE_SE=${searray[$a]}"
   export alien_CLOSE_SE=${searray[$a]}
   echo "main: alien_CLOSE_SE=$alien_CLOSE_SE"
   # source environment before starting copy procedure 
   proxyrenewal 
     
   array=""
   c=1

   echo "main: reading from file $errorfilesbackup"
   echo "main: checking files"
   if [ -e "$errorfiles" ]; then
    echo "main: $errorfiles existing"
   else
    echo "main: $errorfiles missing"
   fi
   if [ -e "$errorfilesbackup" ]; then
    echo "main: $errorfilesbackup existing"
   else
    echo "main: $errorfilesbackup missing"
   fi


   while read line || [ "$line" != "" ]
    do
     array[c]=$line
     echo "main: array[$c] = $line"
     copyerrorfile
     checkstreams
     c=$(($c+1)) 
     done < $errorfilesbackup

   # wait until old copy processes are terminated
   waitfortermination

  
   errorfilesbackup="$errorfiles-backup$a"
   cp $errorfiles $errorfilesbackup
   echo "main: copying $errorfiles to $errorfilesbackup" 
   echo "main: check if $errorfiles and $errorfilesbackup exist"
   ls $errorfiles
   ls $errorfilesbackup   
   rm $errorfiles
   echo "main: removing old file $errorfiles"
   echo "main: check if $errorfiles and $errorfilesbackup exist"
   ls $errorfiles
   ls $errorfilesbackup

 done

 # copy at least 2 times
 if [ $secount -lt 3 ]; then
   echo "main: one additional copy due to not enough replicas"
   proxyrenewal
   array=""
   c=1
   while read line || [ "$line" != "" ]
    do
     array[c]=$line
     echo "main: array[$c] = $line"
     copyerrorfile
     checkstreams
     c=$(($c+1)) 
     done < $errorfilesbackup
    
   # wait until old copy processes are terminated
   waitfortermination
   
 fi
fi

echo "main: listmode copy skipped"
echo "main: copymode = $copymode"
echo "main: errorfile = $errorfile"

# summary check
echo "main: summary check"

if [ "$config" != "/lustre/alice/alien/$INPUTDIR/" ] ; then
 filearray=($(more $config.out | grep element_count))
 echo "main: DEBUG config = $config"
 echo "main: DEBUG filearray = $filearray"
 filearray_count=${#filearray[*]}
 echo "main: DEBUG filearray_count = $filearray_count"
 filearray_sum=0
 for ((c=2; c<$filearray_count; c=c+3))
  do
  # echo "filearray[$c] = ${filearray[c]}"
  (( filearray_sum= $filearray_sum + ${filearray[c]} ))
 done
 echo "main: files to be copied: $filearray_sum"
 filescopied=($(more $config.out | grep HURRAY -c))
 echo "main: files copied = $filescopied"
 filesskipped=($(more $config.out | grep Skipping -c))
 echo "main: files skipped = $filesskipped"
 copyfailure=($(more $config.out | grep adding -c))
 echo "main: transfer failures = $copyfailure"
 (( sum = $filescopied + $filesskipped + $copyfailure))
 echo "main: sum = $sum"

 # plausibility check
 if [ "$sum" -lt  "$filearray_sum" ] ; then
  echo "warning: files copied + failures + left outs ($sum) are smaller than files in list ($filearray_sum)"
 fi
 
 # check for number of copied files
 (( copysuccess = $filescopied + filesskipped ))
 (( copythreshold = $filearray_sum  - $filearray_sum/10 ))
 echo "copy threshold = $copythreshold must have been copied"
 echo "successful transfers: $copysuccess"
 if [ "$copysuccess" -lt  "$copythreshold" ] ; then
  echo "warning: less files copied ($copysuccess) then specified in threshold ($copythreshold) !!!"
 fi

fi


echo "main: copying $errorfilesbackup back to $errorfiles"
cp $errorfilesbackup $errorfiles

echo "main: removing all wrongly copied files which are still in the list from lustre"
removeerrorfiles

# copy errorfile and config to the processed dir
processeddir="/lustre/alice/alien/$PROCESSEDDIR"
echo "main: moving $errorfiles and $config to $processeddir"
mv $errorfiles $processeddir
chgrp --quiet alidata "$processeddir/$errorfiles"

if [ "$config" != "/lustre/alice/alien/$INPUTDIR/" ] ; then
 mv $config $processeddir 
fi
  
echo "main: done"   

