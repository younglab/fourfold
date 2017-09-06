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


args <- commandArgs(T)

if(length(args) < 2 ) 
  stop ("Not enough arguments")

qnorm <- read.table(paste0(args[1],'/quantile-normalized-samples.txt'),header=T,sep='\t')
qnorm.g <- GRanges(seqnames=as.character(qnorm[,1]),ranges=IRanges(qnorm[,2],width=1),strand='*',as.matrix(qnorm[,-(1:2)]))

region <- convert.string.to.Granges(args[2])

qnorm.g.s <- subsetByOverlaps(qnorm.g,region)

qnorm.s <- as.matrix(mcols(qnorm.g.s))

df <- data.frame(Chromosome=seqnames(qnorm.g.s),Position=start(qnorm.g.s),qnorm.s)

write.table(df,file=args[3],sep='\t',quote=F,row.names = F)

