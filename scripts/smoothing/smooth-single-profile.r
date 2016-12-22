library(matrixStats)
library(GenomicRanges)

write.4c.table <- function(fname,chrs,pos,m) {
  storage.mode(m) <- "double"
  .Call("fast_write",as.character(fname),as.character(chrs),as.integer(pos),m,as.integer(dim(m)))
}

mean.proc <- function(idx,m) {
  
  if(length(idx)>1) {
    p <- m[idx,]
    r <- as.vector(colMeans(m[,-c(1:3)]))
  } else {
    r <- m[idx,-c(1:3)]
  }
  
  r
}

mean.proc.c <- function(fname,chrs,pos,idxl,m) {
  if(ncol(m) > 4 ) {
    m <- m[,-c(1:3)]
  } else {
    m <- matrix(m[,4],ncol=1)
  }
  storage.mode(m) <- "double"
  idxl <- lapply(idxl,as.integer)
  .Call("fourc_smoothing_mean",as.character(fname),as.character(chrs),as.integer(pos),idxl,m,as.integer(dim(m)))
}

linear.proc <- function(idx,m) {
  
  if(length(idx)>1) {
    p <- m[idx,]
    pos <- p[1,1]
    l <- lm.fit(p[,2:3],p[,4]) ## Bx+a model
    
  } else {
    p <- m[idx,]
    pos <- p[1]
    l <- lm.fit(matrix(p[2:3],ncol=2),p[4]) ## Bx+a model
  }

  cf <- coefficients(l)
  if(any(is.na(cf))) cf[is.na(cf)] <- 0
  
  as.vector(cf[2]*pos+cf[1])
}

args <- commandArgs(T)


if(length(args) < 9) {
  stop("failed")
}

library.file <- args[1]
binsize <- as.integer(args[2])
stepsize <- as.integer(args[3])
csize.file <- args[4]
smoothing.mode <- args[5]
measuredsignal.file <- args[6]
bootstrap.file <- args[7]
output.file <- args[8]
output.bootstrap.file <- args[9]


dyn.load(library.file) ### 


proc <- switch(smoothing.mode,
               mean=mean.proc.c)#,
               #linear=linear.proc)

write("debug 0\n",file=stderr())


chrom.sizes <- read.table(csize.file,sep='\t')

signal <- read.table(measuredsignal.file,sep='\t')
bootstrap <- read.table(bootstrap.file,sep='\t')

write("debug 1\n",file=stderr())

s <- GRanges(seqnames=as.character(signal[,1]),ranges=IRanges(signal[,2],width=1),strand='*')#,signal[,3],bootstrap[,-(1:2)])
#d <- cbind(signal[,3],as.matrix(bootstrap[,-c(1:2)]))

### clean up some memory as the next steps will take a lot

g <- unlist(GRangesList(lapply(1:nrow(chrom.sizes),function(i) {
  chr <- as.character(chrom.sizes[i,1])
  msize <- chrom.sizes[i,2]
  
  pos <- seq(1,msize+stepsize,by=stepsize)
  
  GRanges(seqnames=chr,ranges=IRanges(pos,width=binsize),strand='*')
})))

write("debug 2\n",file=stderr())


o <- findOverlaps(g,s)

#m <- cbind(start(g)[queryHits(o)],1,start(s)[subjectHits(o)],d[subjectHits(o),])
ms <- cbind(1,1,start(s),signal[,3])
mb <- cbind(1,1,start(s),as.matrix(bootstrap[,-c(1:2)]))

rm(signal,bootstrap)
gc()


write("debug 3\n",file=stderr())

idx <- split((1:nrow(ms))[subjectHits(o)],queryHits(o))

gs <- g[as.integer(names(idx))]
rm(g)

chrs <- as.character(seqnames(gs))
pos <- as.integer(start(gs))

#r <- sapply(idx,proc,m=m)
#r <- do.call(rbind,lapply(idx,proc,m=m))
proc(output.file,chrs,pos,idx,ms)
proc(output.bootstrap.file,chrs,pos,idx,mb) 


write("debug 4",file=stderr())
#save(r,file=output.file)

