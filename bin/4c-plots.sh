#!/bin/bash

BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/plots
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
INPUTDIR=bootstrap

if [ $# -lt 3 ];
then
  echo "Syntax: 4c-plots.sh <template sample XLSX file> <output dir>"
  exit 1
fi

SAMPLETABLE="$1"
OUTPUTDIR="$2"

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

if [ $? -ne 0 ];
then
  echo "There was an error creating plots, please see error messages"
  exit 1
fi

