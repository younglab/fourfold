#!/bin/bash

BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/plots
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
INPUTDIR=bootstrap

if [ $# -lt 2 ];
then
  echo "Syntax: 4c-plots.sh <template sample XLSX file> <genomic coordinates> <output dir>"
  exit 1
fi

SAMPLETABLE="$1"
COORDINATES="$2"
OUTPUTDIR="$3"

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Cannot find $SAMPLETABLE!"
  exit 1
fi

if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
fi

if [ ! -d "$OUTPUTDIR" ];
then
  echo "$OUTPUTDIR is not a directory!";
  exit 1
fi


# command
$SCRIPTDIR/make-plots.pl $SAMPLETABLE $ORGANISMDATABASE $SCRIPTDIR $COORDINATES $INPUTDIR $OUTPUTDIR

if [ $? -ne 0 ];
then
  echo "There was an error creating plots, please see error messages"
  exit 1
fi

