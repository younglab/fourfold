library(GenomicRanges)
library(matrixStats)

convert.to.granges <- function(fname) {
  temp <- read.table(fname,sep='\t')
  
  GRanges(seqnames=as.character(temp[,1]),ranges=IRanges(temp[,2],width=1),strand='*',signal=temp[,3])
}

build <- function(g,r) {
  v <- rep(0,length(r))
  o <- findOverlaps(r,g)
  
  v[queryHits(o)] <- g$signal[subjectHits(o)]
  
  v
}

args <- commandArgs(T)

if(length(args) < 5) {
  stop ("Need files!")
}

outputfile <- args[1]
pseudocount <- as.numeric(args[2])
args <- args[-(1:2)]

w <- which(args=="SEP")

if(length(w)!=1) {
  stop ("Need to include separator")
}


group1 <- args[1:(w-1)]
group2 <- args[(w+1):length(args)]

g1 <- GRangesList(lapply(group1,convert.to.granges))
g2 <- GRangesList(lapply(group2,convert.to.granges))


r <- reduce(c(unlist(g1),unlist(g2)))

g1.m <- do.call(cbind,lapply(g1,build,r=r))
g2.m <- do.call(cbind,lapply(g2,build,r=r))

g1.mean <- rowMeans(g1.m)
g2.mean <- rowMeans(g2.m)

s.ratio <- log2((g2.mean+pseudocount)/(g1.mean+pseudocount))
#s.bounds <- t(apply(cbind(g1.m,g2.m),1,function(r) {
#  idx <- 1:ncol(g1.m)
#  N <- ncol(g1.m)+ncol(g2.m)
#  m1 <- mean(r[idx])
#  m2 <- mean(r[-idx])
#  s <- sd(r)/sqrt(N)
#  ratio <- (m2+pseudocount)/(m1+pseudocount)
#  x <- qt(.975,N-1)
#  c(log2(ratio-x*s),log2(ratio+x*s))
#}))

#w <- apply(s.bounds,1,function(r) any(!is.finite(r)) || sum(abs(r)>0) < 2 )
w <- apply(cbind(g1.m,g2.m),1,function(r) !all(r>0))

if(any(w)) s.ratio[w] <- 0 ### set anything where the 95% T interval crosses 0 to 0 (too variable to be significant)

mcols(r)[,'ratio'] <- s.ratio


write.table(data.frame(seqnames(r),start(r),r$ratio),outputfile,sep='\t',row.names=F,col.names=F,quote=F)

