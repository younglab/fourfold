#!/bin/bash

BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/quantification

if [ $# -lt 2 ]
then
  echo "Not enough arguments"
  exit 1
fi

TARGETDIR="$1"
REGION="$2"
TARGETOUTPUT="$3"

Rscript $SCRIPTDIR/quantify.r $TARGETDIR $REGION $TARGETOUTPUT

