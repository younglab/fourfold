#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -lt 2 ];
then
  echo "process-4c.sh <sample XSLX table> <bowtie index>"
  exit 1
fi

SAMPLETABLE="$1"
BOWTIEIDX="$2"

### process reads

echo "Processing reads..."

$BASEDIR/process-4c-reads.pl $SAMPLETABLE new

if [ $? -ne 0 ];
then
  echo "Errors in processing the reads, see error messages"
  exit $?
fi


### Align reads
### identify fragments
### map to fragments

echo "Done"
