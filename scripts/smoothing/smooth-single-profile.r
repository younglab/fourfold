library(matrixStats)
library(GenomicRanges)

write.4c.table <- function(fname,chrs,m) {
  storage.mode(m) <- "double"
  .Call("fast_write",as.character(fname),as.character(chrs),m,as.integer(dim(m)))
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

mean.proc.c <- function(idxl,m) {
  m <- m[,-c(1:3)]
  storage.mode(m) <- "double"
  .Call("fourc_smoothing_mean",lapply(idxl,as.integer),m,as.integer(dim(m)))
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


if(length(args) < 8) {
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
d <- cbind(signal[,3],as.matrix(bootstrap[,-c(1:2)]))

### clean up some memory as the next steps will take a lot
rm(signal,bootstrap)
gc()

g <- unlist(GRangesList(lapply(1:nrow(chrom.sizes),function(i) {
  chr <- as.character(chrom.sizes[i,1])
  msize <- chrom.sizes[i,2]
  
  pos <- seq(1,msize+stepsize,by=stepsize)
  
  GRanges(seqnames=chr,ranges=IRanges(pos,width=binsize),strand='*')
})))

write("debug 2\n",file=stderr())


o <- findOverlaps(g,s)

#m <- cbind(start(g)[queryHits(o)],1,start(s)[subjectHits(o)],d[subjectHits(o),])
m <- cbind(1,1,start(s),d)

write("debug 3\n",file=stderr())

idx <- split((1:nrow(m))[subjectHits(o)],queryHits(o))

#r <- sapply(idx,proc,m=m)
#r <- do.call(rbind,lapply(idx,proc,m=m))
r <- proc(idx,m) 

write("debug 4",file=stderr())

#print(r)
#save(r,file='tmp.Rdata')


gs <- g[as.integer(names(idx))]

### clean up some memory as the next step can take a lot of memory
rm(g,m,o,idx)
gc() 

#df <- data.frame(seqnames(gs),start(gs),r)
#df <- cbind(as.character(seqnames(gs)),start(gs),r)

chrs <- as.character(seqnames(gs))
dm <- cbind(start(gs),r)

write.4c.table(output.file,chrs,dm)

#write.table(df,file=output.file,sep='\t',row.names=F,col.names=F,quote=F)
