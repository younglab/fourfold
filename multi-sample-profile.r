library(matrixStats)
library(GenomicRanges)

args <- commandArgs(T)

if(length(args)<8) {
  stop("Need files")
}

output.file <- args[1]
bin.size <- as.integer(args[2])
step.size <- as.integer(args[3])
chrom.file <- args[4]
args <- args[-(1:4)]

chrom.sizes <- read.table(chrom.file,sep='\t')

signal.files <- args[seq(1,length(args),by=2)]
bootstrap.files <- args[seq(2,length(args),by=2)]

N <- length(signal.files)

s <- lapply(signal.files,read.table,sep='\t')
b <- lapply(bootstrap.files,read.table,sep='\t')

gs <- lapply(s,function(df) GRanges(seqnames=as.character(df[,1]),ranges=IRanges(df[,2],width=1),strand='*',signal=df[,3]))

uniq.pos <- reduce(unlist(GRangesList(gs)))

s.v <- lapply(gs,function(g) {
  v <- rep(0,length(uniq.pos))
  
  o <- findOverlaps(uniq.pos,g)
  v[queryHits(o)] <- g$signal[subjectHits(o)]
  v
})

names(s.v) <- paste("sample",1:length(s.v),sep='')

mcols(uniq.pos) <- cbind(mcols(uniq.pos),DataFrame(do.call(cbind,s.v)))

bins <- unlist(GRangesList(lapply(1:nrow(chrom.sizes),function(i) {
  n <- chrom.sizes[i,2]
  
  pos <- seq(0,n+step.size,by=step.size)
  GRanges(seqnames=chrom.sizes[i,1],ranges=IRanges(pos,width=bin.size),strand='*',signal=rep(0,length(pos)))
})))

o <- findOverlaps(bins,uniq.pos)

ms <- sapply(split(as.matrix(mcols(uniq.pos))[subjectHits(o),],queryHits(o)),function(m) {
  if(is.vector(m)) {
    mean(m)
  } else {
    mean(m)
  }
})

idx <- as.integer(names(ms))

df <- data.frame(seqnames(bins)[idx],start(bins)[idx],ms)

write.table(df,file=output.file,sep='\t',row.names=F,col.names=F,quote=F)

