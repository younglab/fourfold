# Fourfold



The details below are meant to be a brief introduction to the pipeline. Expanded instructions will be available at https://github.com/younglab/fourfold/wiki.

## Installation

The package contains a configure script to perform very basic checking of install packages but to also set up the necessary files. (Use autoconf if configure is not present).

To install:

1. Run ./configure --prefix=/path/to/install/directory
2. Run make
3. Run make install

## Template Files

## Running the pipeline

### Launching

The pipeline is executed by running /path/to/install/directory/4c-pipeline.sh <data-sample-template.xslx> <pipeline-template.xslx>. 

If you are running on a LSF-capable cluster, then you can automatically dispatch the pipeline to the cluster via /path/to/install/directory/4c-pipeline.sh -l <data-sample-template.xslx> <pipeline-template.xslx>.

### Troubleshooting

1. The pipeline gets stuck when resuming the analysis

The pipeline tries to resume from the last point it left off in the analysis for all steps smoothing the data and before. This is to speed up the pipeline so it 
doesn't have to repeat a number of time intesive steps. However, if it crashes, the pipeline may get stuck in a bad state. Running the 4c-pipeline.sh frontend
with the "-a" option tells the pipeline to start from the beginning. Additionally, removing all the folders in the output directory would accomplish the same
result.

2. Other general questions/comments

Please e-mail young_computation@wi.mit.edu
