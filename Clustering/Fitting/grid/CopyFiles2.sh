#!/bin/bash

if [[ -z $1 ]]; then
  echo "ERROR: Command line arguments missing. Syntax: CopyFiles.sh [location] [grid software package (optional)]"
  exit 1
fi

# setup the copy for local or grid location
CMD=cp
PREFIX=""
SUFFIX=""
FILEDIR="$(dirname $0)"
FILES="../SelectClusters.C ../CCDBUtils.h ../DataUtils.h ../DigitUtils.h ../PreClusterUtils.h ../TrackUtils.h SelectClusters2.sh"
if [[ $1 == "alien://"* ]]; then
  CMD=alien_cp
  PREFIX="file:"
  SUFFIX="@disk=2"
  FILES+=" merge.C merge.sh validation.sh"
elif [[ ! -d $1 ]]; then
  mkdir $1
fi

# copy the files to the working directory
for file in $FILES; do
  echo "copying $PREFIX$FILEDIR/$file to $1/$(basename $file)"
  $CMD $PREFIX$FILEDIR/$file $1/$(basename $file)$SUFFIX
done

# produce the .jdl files if needed and copy them as well
if [[ $1 == "alien://"* ]]; then

  PACKAGE=VO_ALICE@O2PDPSuite::daily-20250214-0000-1
  if [[ ! -z $2 ]]; then
    PACKAGE=$2
  fi

  WORKINGDIR=${1#"alien://"}
  [[ ! $WORKINGDIR == "/"* ]] && WORKINGDIR=$(alien_home)$WORKINGDIR

  for file in "SelectClusters2.jdl" "merge.jdl" ; do
    sed -e "s#__WORKINGDIR__#$WORKINGDIR#g" -e "s/__PACKAGE__/$PACKAGE/g" $FILEDIR/$file > __$file
    echo "producing and copying __$file to $1/$file"
    $CMD $PREFIX"__"$file $1/$file$SUFFIX
    rm -f __$file
  done

  # also copy the collection if available
  if [[ -f "data.xml" ]]; then
    echo "copying "$PREFIX"data.xml to $1/data.xml"
    $CMD $PREFIX"data.xml" $1/data.xml$SUFFIX
  else
    echo "WARNING: data collection \"data.xml\" not found."
    echo "         produce it and copy it by hand to the working grid directory with the name \"data.xml\""
  fi
fi
