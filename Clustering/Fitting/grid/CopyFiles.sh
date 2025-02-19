#!/bin/bash

if [[ -z $1 ]]; then
  echo "ERROR: Command line arguments missing. Syntax: CopyFiles.sh [location] [grid software package (optional)]"
  exit 1
fi

# setup the copy for local or grid location
CMD=cp
PREFIX=""
FILEDIR="$(dirname $0)"
FILES="merge.C ../SelectClusters.C ../CCDBUtils.h ../DataUtils.h ../DigitUtils.h ../PreClusterUtils.h ../TrackUtils.h SelectClusters.sh"
if [[ $1 == "alien://"* ]]; then
  CMD=alien_cp
  PREFIX="file:"
  FILES+=" merge.sh validation.sh"
elif [[ ! -d $1 ]]; then
  mkdir $1
fi

# copy the files to the working directory
for file in $FILES; do
  echo "copying $PREFIX$FILEDIR/$file to $1/$(basename $file)"
  $CMD $PREFIX$FILEDIR/$file $1/$(basename $file)
done

# produce the .jdl files if needed and copy them as well
if [[ $1 == "alien://"* ]]; then

  PACKAGE=VO_ALICE@O2PDPSuite::daily-20250214-0000-1
  if [[ ! -z $2 ]]; then
    PACKAGE=$2
  fi

  WORKINGDIR=${1#"alien://"}
  [[ ! $WORKINGDIR == "/"* ]] && WORKINGDIR=$(alien_home)$WORKINGDIR

  for file in "SelectClusters.jdl" "merge.jdl" ; do
    sed -e "s#__WORKINGDIR__#$WORKINGDIR#g" -e "s/__PACKAGE__/$PACKAGE/g" $FILEDIR/$file > __$file
    echo "producing and copying __$file to $1/$file"
    $CMD $PREFIX"__"$file $1/$file
    rm -f __$file
  done
fi
