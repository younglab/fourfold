#!/bin/bash

BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/plots
LIBDIR=$BASEDIR/../lib
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
INPUTDIR=bootstrap
SHADING=ci
GROUPONLY=0
YLIMLOW=NA
YLIMHIGH=NA
ENHANCERFILE=NA
PROMOTERFILE=NA


function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  echo "Syntax: 4c-multipe-sample-plots.sh [options] <template sample XLSX file> <multiple sample plot template XSLX file> <genomic coordinates> <output dir>"
  echo "-h help menu"
  echo "-i DIR, --inputdir=DIR set input directory"
  echo "-s TYPE, --shading=TYPE set the shading type, one of confidence interval (ci), standard deviation (sd), or none (na)"
  #echo "--group-only Skip output of individual sample plots"
  echo "--ylim-low=[numeric] Set the lower y-axis range value (must also set --ylim-high)"
  echo "--ylim-high=[numeric] Set the upper y-axis range value (must also set --ylim-low)"
  echo "--add-enhancers=[BED file] Plot enhancers within plot in given BED file"
  #echo "--add-promoters=[BED file] Plot TSS of genes within BED file"
}

TEMP=`getopt -o hi: -l inputdir:,group-only,ylim-low:,ylim-high:,add-enhancers:,add-promoters: -n '4cmultiplots' -- "$@"`
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
	  -i|--inputdir)
	    INPUTDIR="$2"
	    if [ ! -d "$INPUTDIR" ];
	    then
	      helpmenu "Error: $INPUTDIR is not a directory!"
	      exit 1
	    fi
	    shift
	    ;;
	  -s|--shading)
	    SHADING="$2"
	    shift
	    ;;
	  --group-only)
	    GROUPONLY=1
	    ;;
	  --ylim-low)
	    YLIMLOW="$2"
	    shift
	    ;;
	  --ylim-high)
	    YLIMHIGH="$2"
	    shift
	    ;;
	  --add-enhancers)
	    ENHANCERFILE="$2"
	    
	    if [ ! -e "$ENHANCERFILE" ];
	    then
	      helpmenu "Error: provided enhancer file $ENHANCERFILE does not exist!"
	      exit 1
	    fi
	    
	    shift
	    ;;
	   --add-promoters)
	    PROMOTERFILE="$2"
	    
	    if [ ! -e "$PROMOTERFILE" ];
	    then
	      helpmenu "Error: provided enhancer file $PROMOTERFILE does not exist!"
	      exit 1
	    fi
	    
	    shift
	    ;;
	esac
	shift
done

if [ $# -lt 4 ];
then
  helpmenu "Not enough arguments"
  exit 1
fi

SAMPLETABLE="$1"
MULTITABLE="$2"
COORDINATES="$3"
OUTPUTDIR="$4"
#shift
#shift
#shift

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Cannot find $SAMPLETABLE!"
  exit 1
fi

if [ ! -e "$MULTITABLE" ];
then
  echo "Cannot find $MULTITABLE!"
  exit 1
fi


if [ ! -e "$OUTPUTDIR" ];
then
  mkdir $OUTPUTDIR
fi

if [ ! -d "$OUTPUTDIR" ];
then
  echo "$OUTPUTDIR is not a directory!";
  exit 1
fi

if [[ ("$YLIMLOW" == "NA" && "$YLIMHIGH" != "NA") || ("$YLIMLOW" != "NA" && "$YLIMHIGH" == "NA") ]];
then
  helpmenu "Must set both --ylim-low and --ylim-high"
  exit 1
fi


perl -I$LIBDIR $SCRIPTDIR/make-multiple-sample-plots.pl $SAMPLETABLE $MULTITABLE $ORGANISMDATABASE $SCRIPTDIR $COORDINATES $SHADING $INPUTDIR $OUTPUTDIR $YLIMLOW $YLIMHIGH $ENHANCERFILE $PROMOTERFILE

if [ $? -ne 0 ];
then
  echo "There was an error creating plots, please see error messages"
  exit 1
fi

