library(matrixStats)
library(GenomicRanges)

args <- commandArgs(T)


mean.proc <- function(idx,m) {
  mean(m[,3])
}

linear.proc <- function(idx,m) {
  
  pos <- m[idx[1],1]
  l <- lm.fit(cbind(1,matrix(m[idx,2],ncol=1)),m[idx,3]) ## Bx+a model
  
  cf <- coefficients(l)
  
  cf[2]*pos+cf[1]
}

if(length(args) < 6) {
  stop("failed")
}

binsize <- as.integer(args[1])
stepsize <- as.integer(args[2])
csize.file <- args[3]
smoothing.mode <- args[4]
measuredsignal.file <- args[5]
bootstrap.file <- args[6]
output.file <- args[7]

proc <- switch(smoothing.mode,
               mean=mean.proc,
               linear=linear.proc)

write("debug 0\n",file=stderr())


chrom.sizes <- read.table(csize.file,sep='\t')

signal <- read.table(measuredsignal.file,sep='\t')
#bootstrap <- read.table(bootstrap.file,sep='\t')

write("debug 1\n",file=stderr())

s <- GRanges(seqnames=as.character(signal[,1]),ranges=IRanges(signal[,2],width=1),strand='*',signal[,3])#,bootstrap[,-(1:2)])

g <- unlist(GRangesList(lapply(1:nrow(chrom.sizes),function(i) {
  chr <- as.character(chrom.sizes[i,1])
  msize <- chrom.sizes[i,2]
  
  pos <- seq(1,msize+stepsize,by=stepsize)
  
  GRanges(seqnames=chr,ranges=IRanges(pos,width=binsize),strand='*')
})))

write("debug 2\n",file=stderr())


o <- findOverlaps(g,s)

m <- cbind(start(g)[queryHits(o)],start(s)[subjectHits(o)],as.matrix(mcols(s)[subjectHits(o),]))

write("debug 3\n",file=stderr())


#r <- lapply(split(m,queryHits(o)),function(ms) {
#  if(is.vector(ms)) return(as.vector(c(ms[1])))#,quantile(ms[-1],probs=c(.025,.975)))))
#  x <- colMeans(ms)
  
#  c(x[,1])#,as.vector(quantile(x,probs=c(.025,.975))))
#})

r <- sapply(split(1:nrow(m),queryHits(o)),proc,m=m)

write("debug 4",file=stderr())


gs <- g[as.integer(names(r))]
df <- data.frame(seqnames(gs),start(gs),r)

write.table(df,file=output.file,sep='\t',row.names=F,col.names=F,quote=F)
