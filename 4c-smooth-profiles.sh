#!/bin/bash

BASEDIR=$(dirname $0)
ORGANISMDATABASE=$BASEDIR/organism-database.txt
INPUTDIR="bootstrap"

TEMP=`getopt -o hi: -l input-dir: -n '4csmoothing' -- "$@"`
eval set -- "$TEMP"

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
	    ;;
	esac
	shift
done

if [ $# -lt 4 ];
then
  echo "4c-smooth-profiles.sh [options] <XSLX table> <output dir> <bin size> <step size> [sample names]"
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

FILES="all"

if [ $# -gt 0 ];
then
  FILES="$@"
fi

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Cannot find $SAMPLETABLE"
  exit 1
fi

if [ ! -d "$INPUTDIR" ];
then
  echo "Cannot find input directory $INPUTDIR"
  exit 1
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

$BASEDIR/profile-smoothing.pl $SAMPLETABLE $BASEDIR $ORGANISMDATABASE $INPUTDIR $OUTPUTDIR $BINSIZE $STEPSIZE $@

if [ $? -ne 0 ];
then
  echo "Profile smoothing had a problem, see error messages"
  exit 1
fi

