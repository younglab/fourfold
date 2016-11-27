#!/bin/bash

BASEDIR=$(dirname $0)
ORGANISMDATABASE=$BASEDIR/organism-database.txt
PSEUDOCOUNT=0.000001

if [ $# -lt 4 ];
then
  echo "Arguments <template> <input dir> <output dir> <file1> [<file2>] [file3...]"
  exit 1
fi

SAMPLES="$1"
shift
INPUTDIR="$1"
shift
OUTPUTDIR="$1"
shift

if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
fi

if [ ! -d "$OUTPUTDIR" ];
then
  echo "$OUTPUDIR is not a directory!"
  exit 1
fi

$BASEDIR/compare-samples.pl $SAMPLES $ORGANISMDATABASE $BASEDIR $PSEUDOCOUNT $INPUTDIR $OUTPUTDIR $@

if [ $? -ne 0 ];
then
  echo "Comparison failed! see error messages"
  exit 1
fi

