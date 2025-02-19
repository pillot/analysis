#!/bin/bash

if [[ -z $1 ]]; then
  echo "ERROR: Command line arguments missing. Syntax: merge.sh [collection of files to merge]"
  exit 1
fi
if [ ! -f $1 ]; then
  echo "ERROR: input file $1 not found locally"
  exit 1
fi

# extract the list of input dirname/filename from the collection
if [[ "${1##*.}" == "xml" ]]; then
  FILELIST=${1/%xml/txt}
  sed -rn 's/.*turl="([^"]*)".*/\1/p' $1 > $FILELIST
elif [[ "$(file "$1")" =~ ': ASCII text' ]]; then
  FILELIST=$1
else
  echo "invalid input (must be an xml collection of a text file)"
  exit 2
fi

# get the filename to merge
FILE="`head -n1 $FILELIST | xargs -L 1 basename`"
if [[ $(grep -v -c -e "/$FILE$" $FILELIST) != 0 ]]; then
  echo "invalid input (must be a list of dirname/filename with the same filename)"
  exit 2
fi

# do the merging
root -l -b -q 'merge.C+("'$FILELIST'","'$FILE'")'
[[ $? -ne 0 ]] && exit 3

exit 0
