library(affy)
library(limma)
#library(GenomicRanges)


#### 

args <- commandArgs(T)

dfile <- args[1]
issignalfile <- args[2]=='yes'

load(dfile)

if(issignalfile) {
  nsm <- normalizeQuantiles(sm)
  save(nsm,uniq.sites,sample.names,file=dfile)
} else {
  bm <- normalizeQuantiles(bm)
  save(bm,uniq.sites,file=dfile)
}
