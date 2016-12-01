#!/bin/bash

BASEDIR=$(dirname $0)
LIBDIR=$BASEDIR/../lib
SCRIPTDIR=$BASEDIR/../scripts/smoothing
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
INPUTDIR="bootstrap"

TEMP=`getopt -o hi: -l input-dir: -n '4csmoothing' -- "$@"`
eval set -- "$TEMP"

function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  echo "Syntax: 4c-smooth-profiles.sh [options] <XSLX table> <output dir> <bin size> <step size> [sample names]"
  echo "If no sample names (from the Excel table) are given, then all samples within the Excel table are smoothed"
  echo "Options:"
  echo "-h  Help menu"
  echo "-i|--input  Set input directory (default bootstrap)"
}

while [ $# -ge 1 ]; do
	case "$1" in
	  --)
	    shift
	    break
	    ;;
	  -i|--input-dir)
	    INPUTDIR="$2"
	    shift
	    ;;
	  -h)
	    helpmenu
	    exit 0
	    ;;
	esac
	shift
done

if [ $# -lt 4 ];
then
  helpmenu "Error: not enough arguments"
  exit 1
fi

SAMPLETABLE="$1"
shift
OUTPUTDIR="$1"
shift
BINSIZE="$1"
shift
STEPSIZE="$1"
shift

FILES="all"

if [ $# -gt 0 ];
then
  FILES="$@"
fi

if [ ! -e "$SAMPLETABLE" ];
then
  helpmenu "Error: Cannot find $SAMPLETABLE"
  exit 1
fi

if [ ! -d "$INPUTDIR" ];
then
  helpmenu "Error: Cannot find input directory $INPUTDIR"
  exit 1
fi

if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
fi

if [ -e "$OUTPUTDIR" -a ! -d "$OUTPUTDIR" ];
then
  helpmenu "Error: $OUTPUTDIR is not a directory!"
  exit 1
fi

if [[ "$BINSIZE" < 1 || "$STEPSIZE" < 1 ]];
then
  helpmenu "Error: bin size and step size must be non-negative integers"
  exit 1
fi

perl -I$LIBDIR $SCRIPTDIR/profile-smoothing.pl $SAMPLETABLE $BASEDIR $ORGANISMDATABASE $INPUTDIR $OUTPUTDIR $BINSIZE $STEPSIZE $@

if [ $? -ne 0 ];
then
  echo "Profile smoothing had a problem, see error messages"
  exit 1
fi

