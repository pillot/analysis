#!/bin/sh

# script to run dsymutil on all libraries located in $1

for lib in $(ls $1/lib*.dylib $1/lib*.so)
do
 echo $lib
 dsymutil $lib
done
