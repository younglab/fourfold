#!/bin/bash


BASEDIR=$(dirname $0)
SCRIPTDIR=$BASEDIR/../scripts/core
LIBDIR=$BASEDIR/../lib
ORGANISMDATABASE=$BASEDIR/../db/organism-database.txt
ENDAFTERVALIDATION=0
SKIPVALIDATION=0
GEO=no
LSFQUEUE=14
BOWTIEK=1
BOWTIEM=1
BOWTIEN=1
BOOTSTRAPITERATIONS=1000

function cleanproc() {
  rm -rf wigfiles stats bamfiles logs bootstrap *.fq.gz .tmp*
}

function helpmenu() {
  if [ $# -gt 0 ];
  then
    echo "$@"
  fi
  
  printf "Syntax: 4c-read-processing.sh [options] <4C sample XLSX table>\n"
  printf "\tThis command takes the 4C-seq samples within the XSLX spreadsheet and will trim, align, and plot the 4C signal in WIG files\n"
  printf "\nThe following options are supported:\n"
  (printf " %s\t%s\n" "-h" "print this help menu and exit"
   printf " %s\t%s\n" "-r [INT], --resampling=[INT]" "sets the number of bootstrap iterations to generate (default 1000)"
   printf " %s\t%s\n" "-n [INT], --bowtie-n=[INT]" "sets -n parameter for bowtie 1 to INT (INT must be positive)"
   printf " %s\t%s\n" "-m [INT], --bowtie-m=[INT]" "sets -m parameter for bowtie 1 to INT (INT must be positive)"
   printf " %s\t%s\n" "-k [INT], --bowtie-k=[INT]" "sets -k parameter for bowtie 1 to INT (INT must be positive)"
   printf " %s\t%s\n" "--validate-table-only" "checks the 4C sample file for correctness and then exits, skips further processing"
   printf " %s\t%s\n" "--skip-table-validation" "skips table validation step"
   printf " %s\t%s\n" "--lsf-queue=STR" "sets the name of the LSF queue to dispatch alignment jobs in parallel on (default 'normal')"
   printf " %s\t%s\n" "--geo" "generates GEO-ready FASTQ files for upload with corresponding MD5SUM files") | column -t -s $'\t'
  printf "\n"
  
}


TEMP=`getopt -o hn:k:m:v:r: -l validate-table-only,lsf-queue:,bowtie-n:,bowtie-m:,bowtie-k:,bowtie-v:,skip-table-validation,geo,resampling: -n '4cpipeline' -- "$@"`
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
	  -r|--resampling)
	    BOOTSTRAPITERATIONS="$2"
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
	  --geo)
	    GEO=yes
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

  perl -I$LIBDIR $SCRIPTDIR/validate-table.pl $SAMPLETABLE $ORGANISMDATABASE $ENDAFTERVALIDATION

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

echo "Processing reads..."

perl -I$LIBDIR $SCRIPTDIR/process-4c-reads.pl $SAMPLETABLE $ORGANISMDATABASE $BOWTIEN $BOWTIEK $BOWTIEM $GEO

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
  samtools index bamfiles/$PREFIX.sorted.bam; \
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

perl -I$LIBDIR $SCRIPTDIR/map-to-fragments.pl $SAMPLETABLE $BASEDIR $SCRIPTDIR $ORGANISMDATABASE $BOOTSTRAPITERATIONS

if [ $? -ne 0 ];
then
  echo "Errors in mapping the reads to restriction digest fragments, see error messages"
  exit 1
fi

echo "Done"
