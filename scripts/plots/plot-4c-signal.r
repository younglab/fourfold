library(Cairo)
library(GenomicRanges)

convert.string.to.Granges <- function(s) {
  r <- regexec("([[:alnum:]]+):([[:alnum:]]+)-([[:alnum:]]+)",s)
  
  if(r[[1]][1]==-1) {
    stop(paste("Invalid coordinate string",r))
  }
  
  m <- regmatches(s,r)
  v <- m[[1]][-1]
  
  GRanges(seqnames=v[1],ranges=IRanges(as.integer(v[2]),as.integer(v[3])),strand='*')
}

make.plot <- function(g) {
  x <- start(g)
  y <- g$signal
  
  plot(x,y,type='l',ylim=quantile(y,probs=c(0,.95)),xlab="Genomic Position",ylab="RPM")
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
#bt <- read.table(bt.file,sep='\t')

st.g <- GRanges(seqnames=as.character(st[,1]),ranges=IRanges(st[,2],width=1),strand='*',signal=st[,3])

pdf(pdf.file,width=12,height=6)
make.plot(subsetByOverlaps(st.g,region))
dev.off()

CairoPNG(png.file,width=1200,height=600)
make.plot(subsetByOverlaps(st.g,region))
dev.off()
