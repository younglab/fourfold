#!/bin/bash

set -e

BIN4C=/nfs/young_ata/dsday_scratch/4cpipe/bin

mkdir example
cd example

### trims, maps 4C reads
bsub -K -J 4c $BIN4C/4c-pipeline.sh k562-example-4c-data.xlsx &
wait

### quantile normalizes profiles, saves to quantile-norm file, only looks at sample IDs beginning with K562
$BIN4C/4c-normalize-samples.sh k562-example-4c-data.xlsx quantile quantile-norm K562.*

### smoothes quantile normalized profiles (remove --input=quantile-norm to smooth the unnormalized profiles)
$BIN4C/4c-smooth-profiles.sh --input=quantile-norm k562-example-4c-data.xlsx smoothed-profiles 10000 50 K562.*

### make individual replicate and per condition plots of K562 super-enhancers in MYC TAD
$BIN4C/4c-plots.sh --inputdir=smoothed-profiles --ylim-low=0 --ylim-high=800 --add-enhancers=../K562-se.bed k562-example-4c-data.xlsx chr8:130500000-130770818 sample-plots K562.* 

### makes joint plot between two K562 conditions
$BIN4C/4c-multiple-sample-plots.sh --inputdir=smoothed-profiles --ylim-low=0 --ylim-high=800 --add-enhancers=../K562-se.bed k562-example-4c-data.xlsx k562-example-profile-comparison.xlsx chr8:130500000-130770818 multiple-sample-plots
