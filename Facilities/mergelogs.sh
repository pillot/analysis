#!/bin/sh

for file in nEvents.txt nBadEventsAccess.txt nBadEventsBranch.txt nBadEventsRawData.txt nBadEventsReadout.txt nBadEventsCTP.txt nBadEvents.txt nChunksRecover.txt; do

  if [ `find * -maxdepth 1 -name "$file" | wc -l` != 0 ]; then

    cat */$file | awk -v n=0 '{n+=$0;} END {print n}' > $file

  fi

done

for file in error.txt warning.txt triggerPb.txt trackerPb.txt clusterPb.txt accessPb.txt branchPb.txt rawDataPb.txt readoutPb.txt CTPPb.txt eventPb.txt recoverPb.txt; do

  if [ `find * -maxdepth 1 -name "$file" | wc -l` != 0 ]; then

    cat */$file > $file

  fi

done

