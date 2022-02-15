#!/bin/bash

##################################################

validateout=$(dirname "$0")
validatetime=$(date)
error=0
if [ -z "$validateout" ]; then
   validateout="."
fi

cd "$validateout" || exit 1
validateworkdir=$(pwd)

echo "*******************************************************"
echo "* Automatically generated validation script           *"

echo "* Time:    $validatetime "
echo "* Dir:     $validateout"
echo "* Workdir: $validateworkdir"
echo "* ----------------------------------------------------*"
ls -la ./
echo "* ----------------------------------------------------*"

##################################################

if [ ! -f stderr ]; then
   error=1
   echo "Validation error cause: no stderr" >> .alienValidation.trace
   echo "* ########## Job not validated - no stderr  ###"
   echo "* ########## Job not validated - no stderr  ###" >> stderr
   echo "Error = $error"
fi
if [ ! -f stdout ]; then
   error=1
   echo "Validation error cause: no stdout" >> .alienValidation.trace
   echo "* ########## Job not validated - no stdout  ###"
   echo "* ########## Job not validated - no stdout  ###" >> stderr
   echo "Error = $error"
fi

segViol=$(grep -Ei "Segmentation violation" stdout stderr)
segFault=$(grep -Ei "Segmentation fault" stdout stderr)
glibcErr=$(grep -Ei '\*\*\* glibc detected \*\*\*' stdout stderr)
except=$(grep -Ei "Exception caught" stdout stderr)

if [ "$segViol" != "" ]; then
   error=1
   echo "Validation error cause: Segment. violation" >> .alienValidation.trace
   echo "* ########## Job not validated - Segment. violation  ###"
   echo "* ########## Job not validated - Segment. violation  ###" >> stderr
   echo "$segViol"
   echo "Error = $error"
fi
if [ "$segFault" != "" ]; then
   error=1
   echo "Validation error cause: Segment. fault" >> .alienValidation.trace
   echo "* ########## Job not validated - Segment. fault  ###"
   echo "* ########## Job not validated - Segment. fault  ###" >> stderr
   echo "$segFault"
   echo "Error = $error"
fi
if [ "$glibcErr" != "" ]; then
   error=1
   echo "Validation error cause: glibc detected" >> .alienValidation.trace
   echo "* ########## Job not validated - *** glibc detected ***  ###"
   echo "* ########## Job not validated - *** glibc detected ***  ###" >> stderr
   echo "$glibcErr"
   echo "Error = $error"
fi
if [ "$except" != "" ]; then
   error=1
   echo "Validation error cause: Exception caught" >> .alienValidation.trace
   echo "* ########## Job not validated - Exception caught  ###"
   echo "* ########## Job not validated - Exception caught  ###" >> stderr
   echo "$except"
   echo "Error = $error"
fi

##################################################

MISSINGFILES=""

for file in $ALIEN_JDL_FILESTOCHECK; do
   if [ ! -f "$file" ]; then
      if [ ! -z "$MISSINGFILES" ]; then
         MISSINGFILES="$MISSINGFILES,"
      fi
      MISSINGFILES="${MISSINGFILES}${file}"
   fi
done

if [ ! -z "$MISSINGFILES" ]; then
   error=1
   echo "Validation error cause: Required file(s) not found in the output: $MISSINGFILES" >> .alienValidation.trace
   echo "* Error: Required file(s) not found in the output: $MISSINGFILES"
   echo "* Error: Required file(s) not found in the output: $MISSINGFILES" >> stderr
fi

##################################################

if [ $error = 0 ]; then
   echo "* ----------------   Job Validated  ------------------*"
#   echo "* === Logs std* will be deleted === "
#   rm -f std*
fi
echo "* ----------------------------------------------------*"
echo "*******************************************************"
cd - || exit 1
exit $error
