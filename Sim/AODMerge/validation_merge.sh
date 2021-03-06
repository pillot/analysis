#!/bin/bash

##################################################
validateout=`dirname $0`
validatetime=`date`
validated="0";
error=0
if [ -z $validateout ]
then
    validateout="."
fi

cd $validateout;
validateworkdir=`pwd`;

echo "*******************************************************" >> stdout
echo "* Automatically generated validation script           *" >> stdout

echo "* Time:    $validatetime " >> stdout
echo "* Dir:     $validateout" >> stdout
echo "* Workdir: $validateworkdir" >> stdout
echo "* ----------------------------------------------------*" >> stdout
ls -la ./ >> stdout
echo "* ----------------------------------------------------*" >> stdout

##################################################

if [ ! -f stderr ] ; then
   error=1
   echo "* ########## Job not validated - no stderr  ###"  >> stdout
   echo "Error = $error"  >> stdout
fi
parArch=`grep -Ei "Cannot Build the PAR Archive" stderr`
segViol=`grep -Ei "Segmentation violation" stderr`
segFault=`grep -Ei "Segmentation fault" stderr`
fltPoint=`grep -Ei "floating point exception" stderr`
glibcErr=`grep -Ei "*** glibc detected ***" stderr`

if [ "$parArch" != "" ] ; then
   error=1
   echo "* ########## Job not validated - PAR archive not built  ###"  >> stdout
   echo "$parArch"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$segViol" != "" ] ; then
   error=1
   echo "* ########## Job not validated - Segment. violation  ###"  >> stdout
   echo "$segViol"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$segFault" != "" ] ; then
   error=1
   echo "* ########## Job not validated - Segment. fault  ###"  >> stdout
   echo "$segFault"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$fltPoint" != "" ] ; then
   error=1
   echo "* ########## Job not validated - Float. point exception  ###"  >> stdout
   echo "$fltPoint"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$glibcErr" != "" ] ; then
   error=1
   echo "* ########## Job not validated - *** glibc detected ***  ###"  >> stdout
   echo "$glibcErr"  >> stdout
   echo "Error = $error"  >> stdout
fi
# Only check all desired files have been merged properly.
# Skip the validation by the analysis manager after the Terminate
# since all the output files have not been registered at the
# previous stage (like EventStat_temp.root or pyxsec_hists.root)
# and the manager will complain about...
if ! [ -f outputs_valid_merge ] ; then
   error=1
   echo "Output files were not validated by the analysis manager" >> stdout
   echo "Output files were not validated by the analysis manager" >> stderr
fi
if [ $error = 0 ] ; then
   echo "* ----------------   Job Validated  ------------------*" >> stdout
#   echo "* === Logs std* will be deleted === "
#   rm -f std*
fi
echo "* ----------------------------------------------------*"
echo "*******************************************************"
cd -
exit $error
