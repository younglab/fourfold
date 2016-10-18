#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -lt 3 ];
then
  echo "process-4c.sh <sample XSLX table> <bowtie index> <genome FASTA>"
  exit 1
fi

SAMPLETABLE="$1"
BOWTIEIDX="$2"
GENOMEFA="$3"

if [ ! -e "$SAMPLETABLE" ];
then
  echo "Error: $SAMPLETABLE does not exist!"
  exit 1
fi

if [ ! -e "$BOWTIEIDX.1.ebwt" ];
then
  echo "Error: $BOWTIEIDX does not exist!"
  exit 1
fi

if [ ! -e "$GENOMEFA" ];
then
  echo "Error: $GENOMEFA does not exist!"
  exit 1
fi

### process reads

echo "Processing reads..."

$BASEDIR/process-4c-reads.pl $SAMPLETABLE new $BOWTIEIDX

if [ $? -ne 0 ];
then
  echo "Errors in processing the reads, see error messages"
  exit 1
fi

for F in *trimmed.fq
do
  bsub -K -J 4calign "bowtie -n 1 $bowtieidx -p 8 -k 1 -m 1 -S --chunkmbs 256 --best --strata $name.trimmed.fq > $name.sam; \
  gzip $name.trimmed.fq; \
  samtools view -Sb $name.sam > $name.bam; \
  samtools sort -@ 6 -Ttmp $name.bam > $name.sorted.bam";
done

wait

### identify fragments

$BASEDIR/re-fragment-identification.pl $SAMPLETABLE $GENOMEFA

if [ $? -ne 0 ];
then
  echo "Errors in identifying the restriction fragments, see error messages"
  exit 1
fi

### map to fragments

mkdir wigfiles
mkdir stats
$BASEDIR/map-to-fragments.pl $SAMPLETABLE

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
