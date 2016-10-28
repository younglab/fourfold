#!/bin/bash

BASEDIR=$(dirname $0)
ORGANISMDATABASE=$BASEDIR/organism-database.txt

if [ $# -lt 5 ];
then
  echo "4c-smooth-profiles.sh <XSLX table> <output dir> <bin size> <step size> <sample names>"
  exit
fi

SAMPLETABLE="$1"
shift
OUTPUTDIR="$1"
shift
BINSIZE="$1"
shift
STEPSIZE="$1"
shift

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Cannot find $SAMPLETABLE"
  exit
fi

if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
fi

if [ -e "$OUTPUTDIR" -a ! -d "$OUTPUTDIR" ];
then
  echo "$OUTPUTDIR is not a directory!"
  exit 1
fi

$BASEDIR/profile-smoothing.pl $SAMPLETABLE $BASEDIR $ORGANISMDATABASE $OUTPUTDIR $BINSIZE $STEPSIZE $@

if [ $? -ne 0 ];
then
  echo "Profile smoothing had a problem, see error messages"
  exit 1
fi

