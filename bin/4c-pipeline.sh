#!/bin/bash

BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/pipeline
LIBDIR=$BASEDIR/../lib
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
RUNALL=0

function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  printf "Syntax: 4c-pipeline.sh [options] <4C pipeline XLSX file>\n"
  printf "\n"
  printf "\nThe following options are supported:\n"
  (printf " %s\t%s\n" "-h" "print this help menu and exits"
   printf " %s\t%s\n" "-a, --run-all" "ignores whether any files already exist, starts the pipeline from the beginning") | column -t -s $'\t'
  printf "\n"
  
}

TEMP=`getopt -o ha -l run-all -n '4cpipeline' -- "$@"`
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
	esac
	shift
done

if [ $# -lt 1 ];
then
  helpmenu "Error: not enough arguments"
  exit 1
fi

PIPELINEFILE="$1"

if [ ! -e "$PIPELINEFILE" ];
then
  helpmenu "Error: $PIPELINEFILE does not exist!"
  exit 1
fi

$SCRIPTDIR/managed-pipeline.pl $BASEDIR $RUNALL $PIPELINEFILE

if [ $? -ne 0 ];
then
  echo "There was an error in running the pipeline, please see error messages"
  exit 1
fi

echo "4C pipeline is finished"
