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
VERTLINES=NA
CIALPHA=0.95


function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  printf "Syntax: 4c-plots.sh [options] <template sample XLSX file> <genomic coordinates> <output dir> <files/pattern...>\n"
  printf "\n"
  (printf " %s\t%s\n" "-h" "print help menu and exit"
  printf " %s\t%s\n" "-i DIR, --inputdir=DIR" "sets input directory to read from (default bootstrap)"
  printf " %s\t%s\n" "-s TYPE, --shading=TYPE" "set the shading type, one of confidence interval (ci), standard deviation (sd), or none (na)"
  printf " %s\t%s\n" "--ci-alpha=[PERC]" "set the alpha for the CI, needs to be in the range (0,1] (default is .95)"
  printf " %s\t%s\n" "--group-only" "skips the plotting of individual replicates per cell type and condition"
  printf " %s\t%s\n" "--ylim-low=[numeric]" "set the value of the bottom of the y-axis range (must also set --ylim-high)"
  printf " %s\t%s\n" "--ylim-high=[numeric]" "set the value of the top of the y-axis range (must also set --ylim-low)"
  printf " %s\t%s\n" "--add-enhancers=[BED file]" "plot any enhancers in the BED file are in the range of the genomic coordinates"
  printf " %s\t%s\n" "--add-vertical-line=[POS1,POS2,...]" "plot black vertical line at positions POS1[,POS2,...], must be comma-delimited") | column -t -s $'\t'
  #echo "--add-promoters=[BED file] Plot TSS of genes within BED file"
  
  printf "\n"
}

TEMP=`getopt -o hi:s: -l inputdir:,group-only,ylim-low:,ylim-high:,add-enhancers:,add-promoters:,shading:,add-vertical-line:,ci-alpha: -n '4cplots' -- "$@"`
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
	  --add-vertical-line)
	    VERTLINES="$2"
	    shift
	    ;;
	  --ci-alpha)
	    CIALPHA="$2"
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
COORDINATES="$2"
OUTPUTDIR="$3"
shift
shift
shift

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Cannot find $SAMPLETABLE!"
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


perl -I$LIBDIR $SCRIPTDIR/make-plots.pl $SAMPLETABLE $ORGANISMDATABASE $SCRIPTDIR $COORDINATES $SHADING $INPUTDIR $OUTPUTDIR $GROUPONLY $YLIMLOW $YLIMHIGH $ENHANCERFILE $PROMOTERFILE $VERTLINES $CIALPHA $@

if [ $? -ne 0 ];
then
  echo "There was an error creating plots, please see error messages"
  exit 1
fi

