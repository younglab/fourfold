#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -lt 3 ];
then
  echo "4c-smooth-profiles.sh <XSLX table> <output dir> <sample names>"
  exit
fi

SAMPLETABLE="$1"
shift
OUTPUTDIR="$1"
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

$BASEDIR/profile-smoothing.pl $SAMPLETABLE $BASEDIR $OUTPUTDIR $@

if [ $? -ne 0 ];
then
  echo "Profile smoothing had a problem, see error messages"
  exit 1
fi

