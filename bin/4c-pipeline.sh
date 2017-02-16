#!/bin/bash

BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/pipeline
LIBDIR=$BASEDIR/../lib
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
RUNALL=0
OUTPUTDIR=.

function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  printf "Syntax: 4c-pipeline.sh [options] <4C data template XLSX file> <4C pipeline XLSX file>\n"
  printf "\n"
  printf "\nThe following options are supported:\n"
  (printf " %s\t%s\n" "-h" "print this help menu and exits"
   printf " %s\t%s\n" "-a, --run-all" "ignores whether any files already exist, starts the pipeline from the beginning"
   printf " %s\t%s\n" "-o, --output" "sets the output directory (default to generate output in current directory)") | column -t -s $'\t'
  printf "\n"
  
}

TEMP=`getopt -o hao: -l run-all,output: -n '4cpipeline' -- "$@"`
eval set -- "$TEMP"

while [ $# -ge 1 ]; do
	case "$1" in
	  --)
	    shift
	    break
	    ;;
	  -h)
	    helpmenu
	    exit 0
	    ;;
	  -a|--run-all)
	    RUNALL=1
	    ;;
	  -o|--output)
	    OUTPUTDIR="$2"
	    shift
	    ;;
	esac
	shift
done

if [ $# -lt 2 ];
then
  helpmenu "Error: not enough arguments"
  exit 1
fi

DATAFILE="$1"
PIPELINEFILE="$2"

if [ ! -e "$DATAFILE" ];
then
  helpmenu "Error: $DATAFILE does not exist!"
  exit 1
fi

if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
  if [ $? -ne 0 ];
  then
    helpmenu "Error: Failed to create $OUTPUTDIR"
    exit 1
  fi
fi

if [ ! -d "$OUTPUTDIR" ];
then
  helpmenu "Error: $OUTPUTDIR is not a directory!"
  exit 1
fi
  
if [ ! -e "$PIPELINEFILE" ];
then
  helpmenu "Error: $PIPELINEFILE does not exist!"
  exit 1
fi

$SCRIPTDIR/managed-pipeline.pl $BASEDIR $RUNALL $DATAFILE $OUTPUTDIR $PIPELINEFILE

if [ $? -ne 0 ];
then
  echo "There was an error in running the pipeline, please see error messages"
  exit 1
fi

echo "4C pipeline is finished"
