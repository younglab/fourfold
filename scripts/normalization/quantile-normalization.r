library(affy)
library(GenomicRanges)
library(limma)

args <- commandArgs(T)

outdir <- args[1]
output <- paste(outdir,"quantile-normalized-samples.txt",sep='/')
args <- args[-1]

if(length(args)<2) {
  stop("need samples")
}

sample.names <- args[seq(1,length(args),2)]
files <- args[seq(2,length(args),2)]
  
df <- lapply(files,read.table,sep='\t')

g <- lapply(df,function(d) {
  GRanges(seqnames=as.character(d[,1]),ranges=IRanges(d[,2],width=1),strand='*',signal=d[,3])
})

uniq.sites <- reduce(unlist(GRangesList(g)))

sm <- do.call(cbind,lapply(g,function(x) {
  v <- rep(0,length(uniq.sites))
  
  o <- findOverlaps(uniq.sites,x)
  
  v[queryHits(o)] <- x$signal[subjectHits(o)]
  v
}))

nsm <- normalizeQuantiles(sm)

out <- data.frame(seqnames(uniq.sites),start(uniq.sites),nsm)

colnames(out) <- c("chr","pos",sample.names)
rownames(out) <- NULL

write.table(out,file=output,col.names=T,row.names=F,quote=F,sep='\t')

lapply(colnames(out)[-(1:2)],function(n) {
  df <- dfrpm <- out[,c("chr","pos",n)]
  dfrpm[,3] <- df[,3]/sum(df[,3])*1e6
  fcounts <- paste(outdir,paste(n,"filtered.counts.txt",sep='.'),sep='/')
  frpm <- paste(outdir,paste(n,"filtered.rpm.txt",sep='.'),sep='/')
  
  write.table(df,file=fcounts,col.names=F,row.names=F,quote=F,sep='\t')
  write.table(dfrpm,file=frpm,col.names=F,row.names=F,quote=F,sep='\t')
})
