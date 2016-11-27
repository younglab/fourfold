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
outputsmoothedfile <- args[2]
csize.file <- args[3]
binsize <- as.integer(args[4])
stepsize <- as.integer(args[5])
pseudocount <- as.numeric(args[6])
args <- args[-(1:6)]

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

rred <- r[r$ratio!=0 & !is.na(r$ratio)]

chrom.sizes <- read.table(csize.file,sep='\t')

#signal <- read.table(measuredsignal.file,sep='\t')
#bootstrap <- read.table(bootstrap.file,sep='\t')


#s <- GRanges(seqnames=as.character(signal[,1]),ranges=IRanges(signal[,2],width=1),strand='*',signal[,3])#,bootstrap[,-(1:2)])

s <- rred

g <- unlist(GRangesList(lapply(1:nrow(chrom.sizes),function(i) {
  chr <- as.character(chrom.sizes[i,1])
  msize <- chrom.sizes[i,2]
  
  pos <- seq(1,msize+stepsize,by=stepsize)
  
  GRanges(seqnames=chr,ranges=IRanges(pos,width=binsize),strand='*')
})))

o <- findOverlaps(g,s)

m <- as.matrix(mcols(s)[subjectHits(o),])

#r <- lapply(split(m,queryHits(o)),function(ms) {
#  if(is.vector(ms)) return(as.vector(c(ms[1])))#,quantile(ms[-1],probs=c(.025,.975)))))
#  x <- colMeans(ms)

#  c(x[,1])#,as.vector(quantile(x,probs=c(.025,.975))))
#})

r <- sapply(split(m,queryHits(o)),mean)

gs <- g[as.integer(names(r))]
df <- data.frame(seqnames(gs),start(gs),r)#,do.call(rbind,r))

write.table(df,file=outputsmoothedfile,sep='\t',row.names=F,col.names=F,quote=F)

