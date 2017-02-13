#!/bin/bash

BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/normalization
LIBDIR=$BASEDIR/../lib
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt

function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  echo "4c-normalize-samples.sh <sample table> <normalization type> <output dir> <samples...>"
  echo "Valid normalization types are: quantile, quantile-cis"
}


if [ $# -lt 4 ];
then
  helpmenu "Error: not enough arguments"
  exit 1
fi

SAMPLETABLE="$1"
NORMTYPE="$2"
OUTPUTDIR="$3"
shift
shift
shift

if [ ! -e "$SAMPLETABLE" ];
then
  helpmenu "Error: cannot find $SAMPLETABLE!"
  exit 1
fi

if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
fi

if [ ! -d "$OUTPUTDIR" ];
then
  helpmenu "Error: $OUTPUTDIR is not a directory"
  exit 1
fi

echo "Running normalization process: $NORMTYPE"

case $NORMTYPE in
  quantile)
    perl -I$LIBDIR $SCRIPTDIR/process-quantile-norm.pl $SAMPLETABLE $SCRIPTDIR $ORGANISMDATABASE $OUTPUTDIR 0 $@
    ;;
  quantile-cis)
    perl -I$LIBDIR $SCRIPTDIR/process-quantile-norm.pl $SAMPLETABLE $SCRIPTDIR $ORGANISMDATABASE $OUTPUTDIR 1 $@
    ;;
    
  *)
    echo "Error: Unknown normalization type"
    exit 1
    ;;
esac

if [ $? -ne 0 ];
then
  echo "Error in normalizing samples, please see error messages"
  exit 1
fi

