#!/bin/sh

# if $3 != "revert": print files and events which are in $1 but not in $2
# if $3 == "revert": print files and events which are in both $1 and $2

if [[ $# < 2 ]]; then
  echo "give the 2 lists of files@events to compare"
  exit 3
fi

nlines=`cat $1 | wc -l`
iline=0
nevtot=0
nevsame=0
nevdiff=0

for line in `cat $1`; do

  let iline++
  let frac=100*iline/nlines
  printf "processing file... %d%%\r" $frac

  file=`echo $line | cut -d '@' -f 1`
  events=`echo $line | grep '@' | cut -d '@' -f 2 | sed 's/,/\ /g'`
  let nevtot="nevtot+`echo $events | wc -w`"

  line2=`grep "$file" $2`

  if [[ $line2 ]]; then

    events2=`echo $line2 | grep '@' | cut -d '@' -f 2`

    if [[ $events && $events2 ]]; then

      events2=",$events2,"
      sameevents=""
      missingevents=""

      for event in $events; do

        if [[ `echo $events2 | grep -v ",$event,"` ]]; then
          missingevents="$missingevents,$event"
          let nevdiff++
        else
          sameevents="$sameevents,$event"
          let nevsame++
        fi

      done

      if [[ $3 != "revert" ]]; then
        if [[ $missingevents ]]; then
          echo "$file@${missingevents/,/}"
        fi
      else
        if [[ $sameevents ]]; then
          echo "$file@${sameevents/,/}"
        fi
      fi

    else

      if [[ $3 == "revert" ]]; then
        echo "$line"
      fi

    fi

  else

    let nevdiff="nevdiff+`echo $events | wc -w`"

    if [[ $3 != "revert" ]]; then
      echo "$line"
    fi

  fi

done

echo "event stat: diff=$nevdiff, same=$nevsame, total=$nevtot"
