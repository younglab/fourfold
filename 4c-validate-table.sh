#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -lt 1 ];
then
  echo "Syntax: 4c-validate-table.sh <sample XSLX table>"
  exit 1
fi

SAMPLETABLE="$1"

$BASEDIR/4c-pipeline.sh --validate-table-only $SAMPLETABLE
