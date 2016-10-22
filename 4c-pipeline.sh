#!/bin/bash


BASEDIR=$(dirname $0)
ORGANISMDATABASE=$BASEDIR/organism-database.txt
ENDAFTERVALIDATION=0
LSFQUEUE=normal


TEMP=`getopt -o h -l validate-table-only,lsf-queue: -n '4cpipeline' -- "$@"`
eval set -- "$TEMP"

while [ $# -ge 1 ]; do
	case "$1" in
	  --)
	    shift
	    break
	    ;;
	  -h)
	    ;;
	  --validate-table-only)
	    ENDAFTERVALIDATION=1
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
  echo "process-4c.sh <sample XSLX table>"
  exit 1
fi

SAMPLETABLE="$1"
#BOWTIEIDX="$2"
#GENOMEFA="$3"
MINFRAGMENTLENGTH=20

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Error: $SAMPLETABLE does not exist!"
  exit 1
fi

if [ $ENDAFTERVALIDATION -eq 0 ];
then
  mkdir wigfiles
  mkdir stats
  mkdir bamfiles
  mkdir logs
fi

echo "Validating sample table..."

$BASEDIR/validate-table.pl $SAMPLETABLE $ORGANISMDATABASE

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

### process reads

echo "4C analysis launced in $PWD by $USER" | mail -s "[4C] Analysis Pipeline Started" dsday@wi.mit.edu

echo "Processing reads..."

$BASEDIR/process-4c-reads.pl $SAMPLETABLE $ORGANISMDATABASE

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

$BASEDIR/add-mapping-stats.pl $SAMPLETABLE

### identify fragments

echo "Identifying fragments..."

$BASEDIR/re-fragment-identification.pl $SAMPLETABLE $ORGANISMDATABASE $MINFRAGMENTLENGTH

if [ $? -ne 0 ];
then
  echo "Errors in identifying the restriction fragments, see error messages"
  exit 1
fi

### map to fragments

echo "Mapping to fragments..."

$BASEDIR/map-to-fragments.pl $SAMPLETABLE $BASEDIR

if [ $? -ne 0 ];
then
  echo "Errors in mapping the reads to restriction digest fragments, see error messages"
  exit 1
fi

### smooth profiles

### comparisons

### generate reports

#mkdir reports

echo "4C analysis finished in $PWD by $USER" | mail -s "[4C] Analysis Pipeline Finished" dsday@wi.mit.edu


echo "Done"
