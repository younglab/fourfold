library(GenomicRanges)

read.signal.file <- function(fname) {
  d <- read.table(fname,sep='\t')
  
  GRanges(seqnames=as.character(d[,1]),ranges=IRanges(d[,2],width=1),strand='*',signal=Rle(d[,3]))
}

read.bootstrap.file <- function(fname) {
  d <- read.table(fname,sep='\t')
  
  GRanges(seqnames=as.character(d[,1]),ranges=IRanges(d[,2],width=1),strand='*',as.matrix(d[,-(1:2)]))
}

signal.covert.to.matrix <- function(gsl,uniq.sites,sample.names,output.file) {
  sm <- do.call(cbind,lapply(gsl,function(x) {
    v <- rep(0,length(uniq.sites))
    
    o <- findOverlaps(uniq.sites,x)
    
    v[queryHits(o)] <- as.vector(x$signal)[subjectHits(o)]
    v
  }))

  save(sm,uniq.sites,sample.names,file=output.file)
  output.file
}

background.covert.to.matrix <- function(gsl,uniq.sites,output.prefix) {
  bsm <- lapply(gsl,function(x) {
    v <- matrix(0,nrow=length(uniq.sites),ncol=ncol(mcols(x)))
    
    o <- findOverlaps(uniq.sites,x)
    
    v[queryHits(o),] <- as.matrix(mcols(x))[subjectHits(o),]
    v
  })
  
  bootstrap.idxes <- 1:ncol(bsm[[1]])
  
  files <- paste(output.prefix,"-",bootstrap.idxes,".Rdata",sep='')

  lapply(bootstrap.idxes,function(i) {
    bm <- do.call(cbind,lapply(bsm,function(m) m[,i]))
    
    save(bm,uniq.sites,file=files[i])
  })
  
  files
}


#####


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

uniq.sites <- reduce(unlist(GRangesList(g)))

if(onlycis) {
  uniq.sites <- uniq.sites[seqnames(uniq.sites)==cischrom]
}

signal.output.file <- paste(outdir,".signal-matrix-tmp.Rdata",sep='/')
background.output.prefix <- paste(outdir,".background-matrix-tmp",sep='/')

signal.covert.to.matrix(g,uniq.sites,sample.names,signal.output.file)
g <- NULL
rm(g)
gc()

background.covert.to.matrix(bg,uniq.sites,background.output.prefix)
