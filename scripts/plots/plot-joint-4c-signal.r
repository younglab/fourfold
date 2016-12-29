library(Cairo)
library(GenomicRanges)
library(matrixStats)

convert.string.to.Granges <- function(s) {
  r <- regexec("([[:alnum:]]+):([[:alnum:]]+)-([[:alnum:]]+)",s)
  
  if(r[[1]][1]==-1) {
    stop(paste("Invalid coordinate string",r))
  }
  
  m <- regmatches(s,r)
  v <- m[[1]][-1]
  
  GRanges(seqnames=v[1],ranges=IRanges(as.integer(v[2]),as.integer(v[3])),strand='*')
}

conf.int <- function(x,y,m,isstatsfile) {
  if( isstatsfile) 
    q <- m[,4:5]
  else
    q <- rowQuantiles(m,probs=c(.025,.975))

  polygon(c(x,rev(x)),c(q[,1],rev(y)),col=rgb(1,0,0,.5))
  polygon(c(x,rev(x)),c(q[,2],rev(y)),col=rgb(1,0,0,.5))
  #polygon(x,q[,1],col=rgb(1,0,0,.5))
  #polygon(x,q[,2],col=rgb(1,0,0,.5))
}

int.sd <- function(x,y,m,isstatsfile) {
  if( isstatsfile ) 
    q <- sqrt(m[,2])
  else
    q <- rowSds(m)
  
  polygon(c(x,rev(x)),c(y+q,rev(y)),col=rgb(1,0,0,.5))
  polygon(c(x,rev(x)),c(y-q,rev(y)),col=rgb(1,0,0,.5))
}

polygon.ignore <- function(x,y,m,isstatsfile) {
  
}

#make.plot <- function(g) {
make.plot <- function(x,y,bm,isstatsfile) {
  #x <- start(g)
  #y <- g$signal
  
  plot(x,y,type='l',ylim=quantile(y,probs=c(0,.95)),xlab="Genomic Position",ylab="RPM")
  
  f <- switch(c.type,
              ci=conf.int,
              sd=int.sd,
              na=polygon.ignore)
  
  #f(x,y,as.matrix(mcols(g)[,-1]))
  f(x,y,bm,isstatsfile)
}

args <- commandArgs(T)

if(length(args)<7) {
  stop("Arguments: <is stats file> <type> <region> <output PDF> <output PNG> [<signal table> <bootstrap table> reps] ")
}


isstatsfile <- as.integer(args[1]) != 0
c.type <- args[2]
region <- convert.string.to.Granges(args[3])
pdf.file <- args[4]
png.file <- args[5]
args <- args[-(1:5)]

st.files <- args[seq(1,length(args),2)]
bt.files <- args[seq(2,length(args),2)]

if(!isstatsfile) {
  stop("Currently only works on data sets with the stats bootstrap file!")
}


st <- lapply(st.files,read.table,sep='\t')
bt <- lapply(bt.files,read.table,sep='\t')

uniq.pos <- reduce(unlist(GRangesList(lapply(st,function(df) GRanges(seqnames=as.character(df[,1]),ranges=IRanges(df[,2],width=1),strand='*')))))


#pos <- GRanges(seqnames=as.character(st[[1]][,1]),ranges=IRanges(st[[1]][,2],width=1),strand='*')
o <- findOverlaps(uniq.pos,region)

st.v <- rowMeans(do.call(cbind,lapply(st,function(df) {
  g <- GRanges(seqnames=as.character(df[,1]),ranges=IRanges(df[,2]),strand='*')
  
  x <- findOverlaps(uniq.pos,g)
  
  v <- rep(0,length(uniq.pos))
  
  v[queryHits(x)] <- df[,3][subjectHits(o)]
  
  v
})))[queryHits(o)]

bt.m <- do.call(cbind,lapply(3:7,function(i) rowMeans(do.call(cbind,lapply(st,function(df) {
  g <- GRanges(seqnames=as.character(df[,1]),ranges=IRanges(df[,2]),strand='*')
  
  x <- findOverlaps(uniq.pos,g)
  
  v <- rep(0,length(uniq.pos))
  
  v[queryHits(x)] <- df[,i][subjectHits(o)]
  
  v
})))))[queryHits(o),]

pdf(pdf.file,width=12,height=6)
make.plot(start(pos)[queryHits(o)],st.v,bt.m,isstatsfile)
dev.off()

CairoPNG(png.file,width=1200,height=600)
make.plot(start(pos)[queryHits(o)],st.v,bt.m,isstatsfile)
dev.off()
