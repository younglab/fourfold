#!/bin/bash

BASEDIR=$(dirname $0)

function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  echo "Syntax: 4c-validate-table.sh <sample XSLX table>"

}

if [ $# -lt 1 ];
then
  helpmenu "Error: not enough arguments"
  exit 1
fi

SAMPLETABLE="$1"

if [ ! -e "$SAMPLETABLE" ];
then
  helpmenu "Error: $SAMPLETABLE does not exist!"
  exit 1
fi

$BASEDIR/4c-read-processing.sh --validate-table-only $SAMPLETABLE
