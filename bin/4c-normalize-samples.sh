#!/bin/bash

BASEDIR=$(dirname $0)
ORGANISMDATABASE=$BASEDIR/organism-database.txt


if [ $# -lt 4 ];
then
  echo "4c-normalize-samples.sh <sample table> <normalization type> <output dir> <samples...>"
  echo "Valid normalization types are: quantile, "
  exit 1
fi

SAMPLETABLE="$1"
NORMTYPE="$2"
OUTPUTDIR="$3"
shift
shift
shift

if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
fi

if [ ! -d "$OUTPUTDIR" ];
then
  echo "Error: $OUTPUTDIR is not a directory"
  exit 1
fi

case $NORMTYPE in
  quantile)
    $BASEDIR/process-quantile-norm.pl $SAMPLETABLE $BASEDIR $ORGANISMDATABASE $OUTPUTDIR $@
    ;;
  
  *)
    echo "Unknown normalization type"
    exit 1
    ;;
esac

