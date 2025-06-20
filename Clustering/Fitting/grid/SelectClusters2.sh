#!/bin/bash

if [[ -z $2 ]]; then
  echo "ERROR: Command line arguments missing. Syntax: SelectClusters.sh [collection of one single file to process] [run number]"
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

# this version works only with a single input
if [[ $(cat $FILELIST | wc -l) -ne 1 ]]; then
  echo "the collection must contain one single input"
  exit 2
fi

# get the location of the input files
DIR="`head -n1 $FILELIST | xargs -L 1 dirname`"

# do the merging of input files
for file in "mchtracks.root" "mchclusters.root" "muontracks.root"; do
  echo "copying $DIR/$file locally"
  alien_cp $DIR/$file file:
done

# run the cluster selection
SELECTCLUSTERS='SelectClusters.C+('$2',"mchclusters.root","mchtracks.root","muontracks.root",true,false,"clusters.root")'
echo $'\n'Processing $SELECTCLUSTERS...
root -l -b <<EOF
gSystem->Load("libO2MCHMappingImpl4");
gSystem->Load("libO2MCHTracking");
.x $SELECTCLUSTERS
EOF
[[ $? -ne 0 ]] && exit 3

exit 0
