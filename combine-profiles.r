library(GenomicRanges)

args <- commandArgs(T)

if(length(args)<4) {
  stop("Need at least 4 arguments: sample type, condition, output directory, file 1, file 2, [file 3...]")
}

sampletype <- args[1]
condition <- args[2]
outdir <- args[3]
files <- args[-c(1:3)]

gl <- GRangesList(lapply(files,function(f) {
  temp <- read.table(f,sep='\t')
  
  GRanges(seqnames=as.character(temp$V1),ranges=IRanges(temp$V2,width=2),signal=temp$V3)
}))

r <- reduce(unlist(gl))

m <- do.call(cbind,lapply(gl,function(g) {
  v <- rep(0,length(r))
  
  o <- findOverlaps(g,r)
  
  v[subjectHits(o)] <- g$signal[queryHits(o)]
  v
}))

ofile <- paste(outdir,paste(paste(sampletype,condition,sep='-'),"txt",sep='.'),sep='/')

df <- data.frame(chr=seqnames(r),pos=start(r),mean.signal=rowMeans(m))

write.table(df,file=ofile,col.names=F,row.names=F,sep='\t',quote=F)
