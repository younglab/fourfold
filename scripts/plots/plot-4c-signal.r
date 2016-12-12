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

conf.int <- function(x,y,m) {
  q <- rowQuantiles(m,probs=c(.025,.975))

  polygon(c(x,rev(x)),c(q[,1],y),col=rgb(1,0,0,.5))
  polygon(c(x,rev(x)),c(q[,2],y),col=rgb(1,0,0,.5))
}

int.sd <- function(x,y,m) {
  q <- rowSds(m)
  
  polygon(c(x,rev(x)),c(y+q,y),col=rgb(1,0,0,.5))
  polygon(c(x,rev(x)),c(y-q,y),col=rgb(1,0,0,.5))
}

polygon.ignore <- function(x,y,m) {
  
}

make.plot <- function(g) {
  x <- start(g)
  y <- g$signal
  
  plot(x,y,type='l',ylim=quantile(y,probs=c(0,.95)),xlab="Genomic Position",ylab="RPM")
  
  f <- switch(c.type,
              ci=conf.int,
              sd=int.sd,
              na=polygon.ignore)
  
  f(x,y,as.matrix(mcols(g)[,-1]))
}

args <- commandArgs(T)

if(length(args)<5) {
  stop("Arguments: <signal table> <bootstrap table> <type> <region> <output PDF> <output PNG>")
}

st.file <- args[1]
bt.file <- args[2]
c.type <- args[3]
region <- convert.string.to.Granges(args[4])
pdf.file <- args[5]
png.file <- args[6]


st <- read.table(st.file,sep='\t')
bt <- read.table(bt.file,sep='\t')

st.g <- GRanges(seqnames=as.character(st[,1]),ranges=IRanges(st[,2],width=1),strand='*',signal=st[,3],as.matrix(bt[,-c(1:2)]))

pdf(pdf.file,width=12,height=6)
make.plot(subsetByOverlaps(st.g,region))
dev.off()

CairoPNG(png.file,width=1200,height=600)
make.plot(subsetByOverlaps(st.g,region))
dev.off()
