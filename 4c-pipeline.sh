#!/bin/bash


BASEDIR=$(dirname $0)
ORGANISMDATABASE=$BASEDIR/organism-database.txt

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


mkdir wigfiles
mkdir stats
mkdir bamfiles

echo "Validating sample table..."

$BASEDIR/validate-table.pl $SAMPLETABLE $ORGANISMDATABASE

if [ $? -ne 0 ];
then
  echo "Errors in validating the sample table, see error messages"
  exit 1
fi

### process reads

echo "Processing reads..."

$BASEDIR/process-4c-reads.pl $SAMPLETABLE $ORGANISMDATABASE

if [ $? -ne 0 ];
then
  echo "Errors in processing the reads, see error messages"
  exit 1
fi

for F in *trimmed.fq
do
  PREFIX=${F%.trimmed.fq}
  bsub -J 4calign "bash $PREFIX.align.sh; \
  gzip $PREFIX.trimmed.fq; \
  samtools view -Sb bamfiles/$PREFIX.sam > bamfiles/$PREFIX.bam; \
  samtools sort -@ 6 -Ttmp$PREFIX bamfiles/$PREFIX.bam > bamfiles/$PREFIX.sorted.bam; \
  rm bamfiles/$PREFIX.sam bamfiles/$PREFIX.bam;";
done

wait

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

echo "Done"
