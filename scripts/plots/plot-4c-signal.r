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

conf.int <- function(x,y,m,isstatsfile,col,alpha) {
  if( isstatsfile) 
    q <- m[,4:5]
  else
    q <- rowQuantiles(m,probs=c(alpha/2,1-alpha/2))

  polygon(c(x,rev(x)),c(q[,1],rev(y)),col=col,border=NA)
  polygon(c(x,rev(x)),c(q[,2],rev(y)),col=col,border=NA)
}

int.sd <- function(x,y,m,isstatsfile,col,alpha) {
  if( isstatsfile ) 
    q <- sqrt(m[,2])
  else
    q <- rowSds(m)
  
  polygon(c(x,rev(x)),c(y+q,rev(y)),col=col,border=NA)
  polygon(c(x,rev(x)),c(y-q,rev(y)),col=col,border=NA)
}

polygon.ignore <- function(...) { #x,y,m,isstatsfile,col) {
  
}

make.plot <- function(pl,sl,ylimlow,ylimhigh,coord) {
  xval <- unlist(pl)
  xlim <- range(xval)
  
  if(!is.na(ylimlow) && !is.na(ylimhigh)) {
    ylim <- c(ylimlow,ylimhigh)
  } else {
    ylim <- quantile(unlist(sl),probs=c(0,.95))
  }
  
  plot(xval,rep(0,length(xval)),type='n',
       xlim=xlim,ylim=ylim,frame.plot = F,
       ylab="4C-seq signal (RPM)",xlab="Genomic Position (bp)",
       main=coord)
  
  return(ylim)
}
  
draw.sample <- function(x,y,bm,isstatsfile,l.col,s.col,a.trans,alpha) {
  lines(x,y,col=l.col,lwd=2)
  
  f <- switch(c.type,
              ci=conf.int,
              sd=int.sd,
              na=polygon.ignore)
  
  v <- as.vector(col2rgb(s.col))
  col <- rgb(v[1]/255,v[2]/255,v[3]/255,a.trans)
  
  f(x,y,bm,isstatsfile,col,alpha)
}

draw.enhancers.promoters <- function(enhancers,prom,ylim) {
  r <- (ylim[2]-ylim[1])*.05
  
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
    #  segments(start(prom)[i],ylim[1],start(prom)[1],r)
    #  arrows(start(prom)[i],r,start(prom)[i]+d,r)
    #}
  }
}

draw.lines <- function(vertlines) {
  if(vertlines=="NA") return(NULL)
  
  i <- as.integer(unlist(strsplit(vertlines,split=',')))
  
  abline(v=i)
}

args <- commandArgs(T)

if(length(args)<14) {
  stop(paste("Arguments: <region> <type> <is stats file> <ylim low> <ylim high> <enhancer file> <promoter file> <verticle line coordinates> <output PDF> <output PNG> <ci alpha> [<signal table> <bootstrap table> <line color> <shading color> <alpha transparency>]x. Saw ",paste(args)))
}

coord.str <- args[1]
region <- convert.string.to.Granges(args[1])
c.type <- args[2]
isstatsfile <- as.integer(args[3]) != 0
ylimlow <- as.integer(args[4])
ylimhigh <- as.integer(args[5])
enhancerfile <- args[6]
promoterfile <- args[7]
vertlines <- args[8]
pdf.file <- args[9]
png.file <- args[10]
ci.alpha <- as.numeric(args[[1]])

args <- args[-(1:10)]

st.file <- args[seq(1,length(args),5)]
bt.file <- args[seq(2,length(args),5)]
l.colors <- args[seq(3,length(args),5)]
s.colors <- args[seq(4,length(args),5)]
a.trans <- sapply(args[seq(5,length(args),5)],function(ch) as.numeric(ch)/100)


st <- lapply(st.file,read.table,sep='\t')
bt <- if( isstatsfile ) lapply(bt.file,read.table,sep='\t',colClasses = c("factor","integer","numeric","numeric","numeric","numeric","numeric")) else lapply(bt.file,read.table,sep='\t')

enhancers <- NA

if(enhancerfile != "NA") {
  temp <- read.table(enhancerfile,sep='\t')
  e <- GRanges(seqnames=as.character(temp[,1]),ranges=IRanges(temp[,2],temp[,3]),strand='*')
  
  x <- subsetByOverlaps(e,region)
  if(length(x)>0) enhancers <- x
}

prom <- NA

if(promoterfile != "NA") {
  temp <- read.table(promoterfile,sep='\t')
  p <- GRanges(seqnames=as.character(temp[,1]),ranges=IRanges(temp[,2],temp[,3]),strand='*')
  p <- promoters(p,1,1)
  
  x <- subsetByOverlaps(p,region)
  if(length(x)>0) prom <- x
}


pos <- lapply(st,function(m) GRanges(seqnames=as.character(m[,1]),ranges=IRanges(m[,2],width=1),strand='*'))
o <- lapply(pos,function(p) findOverlaps(p,region))
signal <- mapply(function(m,p) m[,3][queryHits(p)],st,o,SIMPLIFY=F)
background <- mapply(function(m,p) as.matrix(m[,-(1:2)])[queryHits(p),],bt,o,SIMPLIFY=F)
allpos <- mapply(function(g,p) start(g)[queryHits(p)],pos,o,SIMPLIFY=F)

pdf(pdf.file,width=12,height=6)
ylim <- make.plot(allpos,signal,ylimlow,ylimhigh,coord.str)
mapply(draw.sample,allpos,signal,background,as.list(rep(isstatsfile,length(allpos))),as.list(l.colors),as.list(s.colors),as.list(a.trans),as.list(rep(alpha,length(allpos))),SIMPLIFY=F)
draw.enhancers.promoters(enhancers,prom,ylim)
draw.lines(vertlines)
dev.off()


CairoPNG(png.file,width=1200,height=600)
ylim <- make.plot(allpos,signal,ylimlow,ylimhigh,coord.str)
mapply(draw.sample,allpos,signal,background,as.list(rep(isstatsfile,length(allpos))),as.list(l.colors),as.list(s.colors),as.list(a.trans),as.list(rep(alpha,length(allpos))),SIMPLIFY=F)
draw.enhancers.promoters(enhancers,prom,ylim)
draw.lines(vertlines)
dev.off()
