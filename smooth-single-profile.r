library(matrixStats)
library(GenomicRanges)

args <- commandArgs(T)

if(length(args) < 6) {
  stop("failed")
}

binsize <- as.integer(args[1])
stepsize <- as.integer(args[2])
csize.file <- args[3]
measuredsignal.file <- args[4]
bootstrap.file <- args[5]
output.file <- args[6]

chrom.sizes <- read.table(csize.file,sep='\t')

signal <- read.table(measuredsignal.file,sep='\t')
bootstrap <- read.table(bootstrap.file,sep='\t')

s <- GRanges(seqnames=as.character(signal[,1]),ranges=IRanges(signal[,2],width=1),strand='*',signal[,3],bootstrap[,-(1:2)])

g <- unlist(GRangesList(lapply(1:nrow(chrom.sizes),function(i) {
  chr <- as.character(chrom.sizes[i,1])
  msize <- chrom.sizes[i,2]
  
  pos <- seq(0,msize+stepsize,by=stepsize)
  
  GRanges(seqnames=chr,ranges=IRanges(pos,width=binsize),strand='*')
})))

o <- findOverlaps(g,s)

m <- as.matrix(mcols(s)[subjectHits(o),])

r <- lapply(split(m,queryHits(o)),function(ms) {
  if(is.vector(ms)) return(as.vector(c(ms[1],quantile(ms[-1],probs=c(.025,.975)))))
  x <- colMeans(ms)
  
  c(x[,1],as.vector(quantile(x,probs=c(.025,.975))))
})

gs <- g[as.integer(names(r))]
df <- data.frame(seqnames(gs),start(gs),end(gs),do.call(rbind,r))

write.table(df,file=output.file,sep='\t',row.names=F,col.names=F,quote=F)
