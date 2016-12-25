library(matrixStats)
library(GenomicRanges)

args <- commandArgs(T)

collapse.ranges <- function(gl,up) {
  v <- rep(0,length(up))
  
  o <- findOverlaps(up,gl[[1]])
  v[queryHits(o)] <- gl[[2]]
  v
}

collapse.ranges.mat <- function(gl,up) {
  v <- matrix(0,length(up),ncol(gl[[2]]))
  
  o <- findOverlaps(up,gl[[1]])
  v[queryHits(o),] <- gl[[2]]
  v
}

mean.proc.c <- function(fname,chrs,pos,idxl,m,stats=0,debug=NULL) {
  if(ncol(m) > 4 ) {
    m <- m[,-c(1:3)]
  } else {
    m <- matrix(m[,4],ncol=1)
  }
  storage.mode(m) <- "double"
  idxl <- lapply(idxl,as.integer)
  
  if(!is.null(debug)) {
    #m <- m[debug,]
    idxl <- idxl[debug]
    #chrs <- chrs[debug]
    #pos <- pos[debug]
  }
  .Call("fourc_smoothing_mean",as.character(fname),as.character(chrs),as.integer(pos),idxl,m,as.integer(dim(m)),as.integer(stats))
}



if(length(args)<10) {
  stop("Need files")
}

library.file <- args[1]
output.file <- args[2]
outputb.file <- args[3]
bin.size <- as.integer(args[4])
step.size <- as.integer(args[5])
chrom.file <- args[6]
args <- args[-(1:6)]

chrom.sizes <- read.table(chrom.file,sep='\t')

signal.files <- args[seq(1,length(args),by=2)]
bootstrap.files <- args[seq(2,length(args),by=2)]

dyn.load(library.file) ### 


N <- length(signal.files)

s <- lapply(signal.files,read.table,sep='\t')
b <- lapply(bootstrap.files,read.table,sep='\t')

ss <- lapply(s,function(df) list(GRanges(seqnames=as.character(df[,1]),ranges=IRanges(df[,2],width=1),strand='*'),matrix(df[,3],ncol=1)) )
bs <- lapply(b,function(df) list(GRanges(seqnames=as.character(df[,1]),ranges=IRanges(df[,2],width=1),strand='*'),as.matrix(df[,-c(1:2)])))

#gs <- lapply(s,function(df) GRanges(seqnames=as.character(df[,1]),ranges=IRanges(df[,2],width=1),strand='*',signal=df[,3]))

uniq.pos <- reduce(unlist(GRangesList(lapply(ss,function(l) l[[1]]))))

#s.v <- lapply(gs,function(g) {
#  v <- rep(0,length(uniq.pos))
#  
#  o <- findOverlaps(uniq.pos,g)
#  v[queryHits(o)] <- g$signal[subjectHits(o)]
#  v
#})

#names(s.v) <- paste("sample",1:length(s.v),sep='')

#mcols(uniq.pos) <- cbind(mcols(uniq.pos),DataFrame(do.call(cbind,s.v)))

s.m <- cbind(1,1,start(uniq.pos),rowMeans(do.call(cbind,lapply(ss,collapse.ranges,up=uniq.pos))))
b.m <- cbind(1,1,start(uniq.pos),do.call(cbind,lapply(bs,collapse.ranges.mat,up=uniq.pos)))

rm(ss,bs)
gc()

bins <- unlist(GRangesList(lapply(1:nrow(chrom.sizes),function(i) {
  n <- chrom.sizes[i,2]
  
  pos <- seq(1,n+step.size,by=step.size)
  GRanges(seqnames=chrom.sizes[i,1],ranges=IRanges(pos,width=bin.size),strand='*')
})))

write("debug 1",file=stderr())


o <- findOverlaps(bins,uniq.pos)

write("debug 2",file=stderr())


idx <- split((1:nrow(s.m))[subjectHits(o)],queryHits(o))

write("debug 3",file=stderr())
print(str(idx))
print(str(s.m))

gs <- uniq.pos[as.integer(names(idx))]

write("debug 4",file=stderr())


#rm(g)

chrs <- as.character(seqnames(gs))
pos <- as.integer(start(gs))

#r <- sapply(idx,proc,m=m)
#r <- do.call(rbind,lapply(idx,proc,m=m))
proc(output.file,chrs,pos,idx,s.m)
proc(outputb.file,chrs,pos,idx,b.m,1) 


write("debug 5",file=stderr())

#ms <- sapply(split(as.matrix(mcols(uniq.pos))[subjectHits(o),],queryHits(o)),mean)
#             #function(m) {
#  if(is.vector(m)) {
#    mean(m)
#  } else {
#    mean(m)
#  }
#})

#idx <- as.integer(names(ms))

#df <- data.frame(seqnames(bins)[idx],start(bins)[idx],ms)

#write.table(df,file=output.file,sep='\t',row.names=F,col.names=F,quote=F)

