#!/bin/bash

MINFRAGMENTLENGTH=20
BASEDIR=$(dirname $0)
ORGANISMDATABASE=$BASEDIR/organism-database.txt

if [ $# -lt 5 ];
then
  echo "Syntax: 4c-run-re-digestion.sh <RE1> <RE1 seq> <RE2> <RE2 seq> <organism>"
  exit 1
fi

RE1="$1"
RE1SEQ="$2"
RE2="$3"
RE2SEQ="$4"
ORGANISM="$5"

echo "Running RE cutting..."

$BASEDIR/standalone-re-cutting.pl $RE1 $RE1SEQ $RE2 $RE2SEQ $ORGANISM

if [ $? -ne 0 ];
then
  echo "Errors in cutting, see error messages"
  exit 1
fi

echo "Finished!"
