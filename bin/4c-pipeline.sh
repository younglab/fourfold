#!/bin/bash


BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/core
LIBDIR=$BASEDIR/../lib
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
ENDAFTERVALIDATION=0
SKIPVALIDATION=0
LSFQUEUE=normal
BOWTIEK=1
BOWTIEM=1
BOWTIEN=1

function cleanproc() {
  rm -rf wigfiles stats bamfiles logs bootstrap *.fq.gz .tmp*
}

function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  echo "Syntax: process-4c.sh <sample XSLX table>"
  echo "-h Help menu (this menu)"
  echo "-n INT, --bowtie-n=INT Sets bowtie 1 -n parameter to INT"
  echo "-m INT, --bowtie-m=INT sets bowtie 1 -m parameter to INT"
  echo "-k INT, --bowtie-k=INT sets bowtie 1 -k parameter to INT"
  echo "--validate-table-only exits after running XSLX table validation"
  echo "--skip-table-validation skips table validation step"
  echo "--lsf-queue=STR sets the LSF queue to dispatch on (default 'normal')"
}


TEMP=`getopt -o hk:m:v: -l validate-table-only,lsf-queue:,bowtie-m:,bowtie-k:,bowtie-v:,skip-table-validation -n '4cpipeline' -- "$@"`
eval set -- "$TEMP"

while [ $# -ge 1 ]; do
	case "$1" in
	  --)
	    shift
	    break
	    ;;
	  -n|--bowtie-n)
	    BOWTIEN="$2"
	    shift
	    ;;
	  -m|--bowtie-m)
	    BOWTIEM="$2"
	    shift
	    ;;
	  -k|--bowtie-k)
	    BOWTIEK="$2"
	    shift
	    ;;
	  -h)
	    helpmenu
	    exit 0
	    ;;
	  --validate-table-only)
	    ENDAFTERVALIDATION=1
	    ;;
	  --skip-table-validation)
	    SKIPVALIDATION=1
	    ;;
	  --lsf-queue)
	    LSFQUEUE="$2"
	    shift
	    ;;
	esac
	shift
done

if [ $# -lt 1 ];
then
  helpmenu "Error: not enough arguments"
  exit 1
fi

SAMPLETABLE="$1"
#BOWTIEIDX="$2"
#GENOMEFA="$3"
MINFRAGMENTLENGTH=20

if [ "$SAMPLETABLE" == "clean" ];
then
  echo "Cleaning directory..."
  cleanproc
  exit 0
fi

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Error: $SAMPLETABLE does not exist!"
  exit 1
fi

if [ "$SKIPVALIDATION" -ne 0 ] && [ "$ENDAFTERVALIDATION" -ne 0 ];
then
  helpmenu "Error: --validate-table-only and --skip-table-validation cannot be both on"
  exit 1
fi

if [ "$ENDAFTERVALIDATION" -eq 0 ];
then
  mkdir wigfiles
  mkdir stats
  mkdir bamfiles
  mkdir logs
  mkdir bootstrap
fi


if [ $SKIPVALIDATION -eq 0 ];
then
  echo "Validating sample table..."

  perl -I$LIBDIR $SCRIPTDIR/validate-table.pl $SAMPLETABLE $ORGANISMDATABASE

  if [ $? -ne 0 ];
  then
    echo "Errors in validating the sample table, see error messages"
    exit 1
  fi

  if [ $ENDAFTERVALIDATION -ne 0 ];
  then
    echo "Sample table is valid"
    exit 0
  fi
else
  echo "Skipping validation step..."
fi

### process reads

echo "4C analysis launced in $PWD by $USER" | mail -s "[4C] Analysis Pipeline Started" dsday@wi.mit.edu

echo "Processing reads..."

perl -I$LIBDIR $SCRIPTDIR/process-4c-reads.pl $SAMPLETABLE $ORGANISMDATABASE $BOWTIEN $BOWTIEK $BOWTIEM

if [ $? -ne 0 ];
then
  echo "Errors in processing the reads, see error messages"
  exit 1
fi

echo "Mapping reads..."

for F in *trimmed.fq
do
  PREFIX=${F%.trimmed.fq}
  bsub -q $LSFQUEUE -o logs/$PREFIX.align.log -K -J 4calign "bash $PREFIX.align.sh; \
  gzip $PREFIX.trimmed.fq; \
  samtools view -Sb bamfiles/$PREFIX.sam > bamfiles/$PREFIX.bam; \
  samtools sort -@ 6 -Ttmp$PREFIX bamfiles/$PREFIX.bam > bamfiles/$PREFIX.sorted.bam; \
  rm bamfiles/$PREFIX.sam bamfiles/$PREFIX.bam;" &
done

wait

perl -I$LIBDIR $SCRIPTDIR/add-mapping-stats.pl $SAMPLETABLE

### identify fragments

echo "Identifying fragments..."

perl -I$LIBDIR $SCRIPTDIR/re-fragment-identification.pl $SAMPLETABLE $ORGANISMDATABASE $MINFRAGMENTLENGTH

if [ $? -ne 0 ];
then
  echo "Errors in identifying the restriction fragments, see error messages"
  exit 1
fi

### map to fragments

echo "Mapping to fragments..."

perl -I$LIBDIR $SCRIPTDIR/map-to-fragments.pl $SAMPLETABLE $BASEDIR $SCRIPTDIR $ORGANISMDATABASE

if [ $? -ne 0 ];
then
  echo "Errors in mapping the reads to restriction digest fragments, see error messages"
  exit 1
fi

### generate reports

#mkdir reports

echo "4C analysis finished in $PWD by $USER" | mail -s "[4C] Analysis Pipeline Finished" dsday@wi.mit.edu


echo "Done"
