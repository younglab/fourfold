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
make.plot <- function(x,y,bm,isstatsfile,ylimlow,ylimhigh,enhancers,prom) {
  
  ylim <- if(!is.na(ylimlow) && !is.na(ylimhigh)) c(ylimlow,ylimhigh) else quantile(y,probs=c(0,.95))
  r <- (ylim[2]-ylim[1])*.05
  
  plot(x,y,type='l',ylim=ylim,xlab="Genomic Position (bp)",ylab="RPM")
  
  if(!is.na(enhancers)) {
    for( i in 1:length(enhancers)) {
      px <- c(start(enhancers)[i],end(enhancers)[i])
      px <- c(px,rev(px))
      
      py <- c(ylim[1],ylim[1],r,r)
      
      polygon(px,py,col='red')
    }
  }
  
  if(!is.na(prom)) {
    #for( i in 1:length(prom)) {
    #  px <- c(start(prom)[i],end(prom)[i])
    #  px <- c(px,rev(px))
    #  
    #  py <- c(ylim[1],ylim[1],r,r)
    #  
    #  polygon(px,py,col='red')
    #}
  }
  
  f <- switch(c.type,
              ci=conf.int,
              sd=int.sd,
              na=polygon.ignore)
  
  f(x,y,bm,isstatsfile)
}

args <- commandArgs(T)

if(length(args)<11) {
  stop("Arguments: <signal table> <bootstrap table> <is stats file> <type> <region> <ylim low> <ylim high> <enhancer file> <promoter file> <output PDF> <output PNG>")
}

st.file <- args[1]
bt.file <- args[2]
isstatsfile <- as.integer(args[3]) != 0
c.type <- args[4]
region <- convert.string.to.Granges(args[5])
ylimlow <- as.integer(args[6])
ylimhigh <- as.integer(args[7])

enhancerfile <- args[8]
promoterfile <- args[9]

pdf.file <- args[10]
png.file <- args[11]


st <- read.table(st.file,sep='\t')
bt <- if( isstatsfile ) read.table(bt.file,sep='\t',colClasses = c("factor","integer","numeric","numeric","numeric","numeric","numeric")) else read.table(bt.file,sep='\t')

pos <- GRanges(seqnames=as.character(st[,1]),ranges=IRanges(st[,2],width=1),strand='*')
o <- findOverlaps(pos,region)

enhancers <- NA

if(enhancerfile != "NA") {
  temp <- read.table(enhancerfile,sep='\t')
  e <- GRanges(seqnames=as.character(temp[,1]),ranges=IRanges(temp[,2],temp[,3]),strand='*')
  
  x <- subsetByOverlaps(e,pos)
  if(length(x)>0) enhancers <- x
}

prom <- NA

if(promoterfile != "NA") {
  temp <- read.table(promoterfile,sep='\t')
  p <- GRanges(seqnames=as.character(temp[,1]),ranges=IRanges(temp[,2],temp[,3]),strand='*')
  
  x <- subsetByOverlaps(p,pos)
  if(length(x)>0) prom <- x
}

st.v <- st[,3][queryHits(o)]
bt.m <- as.matrix(bt[,-(1:2)])[queryHits(o),]

pdf(pdf.file,width=12,height=6)
make.plot(start(pos)[queryHits(o)],st.v,bt.m,isstatsfile,ylimlow,ylimhigh,enhancers,prom)
dev.off()

CairoPNG(png.file,width=1200,height=600)
make.plot(start(pos)[queryHits(o)],st.v,bt.m,isstatsfile,ylimlow,ylimhigh,enhancers,prom)
dev.off()
