#!/bin/bash

function findprimer() {
  echo "Not yet implemented"
  exit 1
}


if [ $# -lt 1 ];
then
  echo "Not enough arguments"
  exit 1
fi

FUNCTION="$1"
shift

case $FUNCTION in
  find-primer)
    TESTFILE="$1"
    
    if [ -z "$TESTFILE" ];
    then
      echo "Need a file to test!"
      exit 1
    fi
    
    if [ ! -e "$TESTFILE" ];
    then
      echo "Cannot find file $TESTFILE"
      exit 1
    fi
  
    findprimer $1
    ;;
  *)
    echo "Unknown utility $FUNCTION"
    exit 1
    ;;
esac

