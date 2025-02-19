#!/bin/bash

if [[ -z $2 ]]; then
  echo "ERROR: Command line arguments missing. Syntax: SelectClusters.sh [collection of files to process] [run number]"
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

# get the filename to be replaced by the input files
FILE0="`head -n1 $FILELIST | xargs -L 1 basename`"
if [[ $(grep -v -c -e "/$FILE0$" $FILELIST) != 0 ]]; then
  echo "invalid input (must be a list of dirname/filename with the same filename)"
  exit 2
fi

# do the merging of input files
for file in "mchtracks.root" "mchclusters.root" "muontracks.root"; do
  sed -rn "s/$FILE0/$file/gp" $FILELIST > "__$FILELIST"
  root -l -b -q 'merge.C+("__'$FILELIST'","'$file'")'
  [[ $? -ne 0 ]] && exit 3
done
rm -f "__$FILELIST"

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
