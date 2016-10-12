#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -lt 3 ];
then
  echo "process-4c.sh <sample XSLX table> <bowtie index> <genome FASTA>"
  exit 1
fi

SAMPLETABLE="$1"
BOWTIEIDX="$2"
GENOMEFA="$3"

### process reads

echo "Processing reads..."

$BASEDIR/process-4c-reads.pl $SAMPLETABLE new $BOWTIEIDX

if [ $? -ne 0 ];
then
  echo "Errors in processing the reads, see error messages"
  exit 1
fi

### identify fragments

$BASEDIR/re-fragment-identification.pl $SAMPLETABLE $GENOMEFA

if [ $? -ne 0 ];
then
  echo "Errors in identifying the restriction fragments, see error messages"
  exit 1
fi

### map to fragments

echo "Done"
