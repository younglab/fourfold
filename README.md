# Fourfold

Fourfold is a 4C-seq processing pipeline. This will handle the processing of 4C-seq reads from FASTQ files on the local disk or the Internet. The output of the pipeline is a set of WIG files, smoothed and unsmoothed, to visualize in a browser like UCSC or IGV as well as PDF plots of selected genomic windows. The PDF plots are a line graph of the signal over the region plus the estimated 95% confidence interval via bootstrap.

The details below are meant to be a brief introduction to the pipeline. Expanded instructions will be available at https://github.com/younglab/fourfold/wiki.

## Installation

Please read INSTALL for installation instructions and dependencies

## Template Files

The necessary template files are in Excel format under the templates/ folder. Each of these needs to be filled out before running the pipeline. Information on what data is needed for each column can be found at https://github.com/younglab/fourfold/wiki.

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

Please file an issue in the tracker on the GitHub page or e-mail young_computation@wi.mit.edu with "[fourfold]" in the subject line.
