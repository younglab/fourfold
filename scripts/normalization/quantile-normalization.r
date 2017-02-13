library(affy)
library(GenomicRanges)
library(limma)

read.signal.file <- function(fname) {
  d <- read.table(fname,sep='\t')
  
  GRanges(seqnames=as.character(d[,1]),ranges=IRanges(d[,2],width=1),strand='*',signal=Rle(d[,3]))
}

read.bootstrap.file <- function(fname) {
  d <- read.table(fname,sep='\t')
  
  GRanges(seqnames=as.character(d[,1]),ranges=IRanges(d[,2],width=1),strand='*',as.matrix(d[,-(1:2)]))
}

signal.normalization <- function(gsl,uniq.sites) {
  sm <- do.call(cbind,lapply(gsl,function(x) {
    v <- rep(0,length(uniq.sites))
    
    o <- findOverlaps(uniq.sites,x)
    
    v[queryHits(o)] <- as.vector(x$signal)[subjectHits(o)]
    v
  }))
  
  normalizeQuantiles(sm)
}

background.normalization <- function(gsl,uniq.sites) {
  bsm <- lapply(gsl,function(x) {
    v <- matrix(0,nrow=length(uniq.sites),ncol=ncol(mcols(x)))
    
    o <- findOverlaps(uniq.sites,x)
    
    v[queryHits(o),] <- as.matrix(mcols(x))[subjectHits(o),]
    v
  })

  
  cat("DEBUG: converted bg into matrices\n")
  
  
  
  do.call(cbind,lapply(1:ncol(bsm[[1]]),function(i) {
    m <- do.call(cbind,lapply(bsm,function(x) x[,i]))
    
    normalizeQuantiles(m)
  }))
  
}

####

args <- commandArgs(T)

print(paste("DEBUG: args length:",length(args)))

outdir <- args[1]
onlycis <- args[2]!="NA"
cischrom <- args[2]

output <- paste(outdir,"quantile-normalized-samples.txt",sep='/')
outputb <- paste(outdir,"quantile-normalized-samples-bootstrap.txt",sep='/')

args <- args[-(1:2)]

if(length(args)<2) {
  stop("need samples")
}

sample.names <- args[seq(1,length(args),3)]
files <- args[seq(2,length(args),3)]
bfiles <- args[seq(3,length(args),3)]

print(paste("DEBUG: outdir:",outdir))
print(paste("DEBUG: onlycis:",onlycis))
print(paste("DEBUG: cischrom:",cischrom))
print(paste("DEBUG: name:",sample.names))
print(paste("DEBUG: files:",files))
print(paste("DEBUG: bfiles:",bfiles))

g <- lapply(files,read.signal.file)
bg <- lapply(bfiles,read.bootstrap.file)

cat("DEBUG: read in files, converted to GRanges\n")

uniq.sites <- reduce(unlist(GRangesList(g)))

if(onlycis) {
  uniq.sites <- uniq.sites[seqnames(uniq.sites)==cischrom]
}

nsm <- signal.normalization(g,uniq.sites)

g <- NULL
rm(g)
gc()

bnsm <- background.normalization(bg,uniq.sites)

bg <- NULL
rm(bg)
gc()

cat("DEBUG: normalized signal\n")


out <- data.frame(seqnames(uniq.sites),start(uniq.sites),nsm)

nsm <- NULL
rm(nsm)
gc()

outb <- data.frame(seqnames(uniq.sites),start(uniq.sites),bnsm)

bnsm <- NULL
rm(bnsm)
gc()

colnames(out) <- c("chr","pos",sample.names)
rownames(out) <- NULL

colnames(outb) <- c("chr","pos",rep(sample.names,ncol(bsm[[1]])))
rownames(outb) <- NULL

cat("DEBUG: converted into table\n")


write.table(out,file=output,col.names=T,row.names=F,quote=F,sep='\t')
write.table(outb,file=outputb,col.names=T,row.names=F,quote=F,sep='\t')

lapply(colnames(out)[-(1:2)],function(n) {
  df <- dfrpm <- out[,c("chr","pos",n)]
  dfb <- dfrpmb <- outb[,colnames(outb) %in% c("chr","pos",n)]
  dfrpm[,3] <- df[,3]/sum(df[,3])*1e6
  dfrpmb[,-c(1:2)] <- apply(as.matrix(dfrpmb[,-c(1:2)]),2,function(r) r/sum(r)*1e6)
  
  fcounts <- paste(outdir,paste(n,"filtered.counts.txt",sep='.'),sep='/')
  frpm <- paste(outdir,paste(n,"filtered.rpm.txt",sep='.'),sep='/')
  
  fcountsb <- paste(outdir,paste(n,"filtered.counts.bootstrap.txt",sep='.'),sep='/')
  frpmb <- paste(outdir,paste(n,"filtered.rpm.bootstrap.txt",sep='.'),sep='/')
  
  
  write.table(df,file=fcounts,col.names=F,row.names=F,quote=F,sep='\t')
  write.table(dfrpm,file=frpm,col.names=F,row.names=F,quote=F,sep='\t')
  write.table(dfb,file=fcountsb,col.names=F,row.names=F,quote=F,sep='\t')
  write.table(dfrpmb,file=frpmb,col.names=F,row.names=F,quote=F,sep='\t')
})
