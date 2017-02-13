library(GenomicRanges)


#######


args <- commandArgs(T)

outdir <- args[1]

output <- paste(outdir,"quantile-normalized-samples.txt",sep='/')
outputb <- paste(outdir,"quantile-normalized-samples-bootstrap.txt",sep='/')


signal.file <- paste(outdir,".signal-matrix-tmp.Rdata",sep='/')
background.files <- dir(path=outdir,pattern="[.]background-matrix",full.names = T,all.files = T)

load(signal.file)


out <- data.frame(seqnames(uniq.sites),start(uniq.sites),nsm)

nsm <- NULL
rm(nsm)
gc()

bnsm <- lapply(background.files,function(f) {
  load(f)
  bm
})

bootstrap.iter <- ncol(bnsm[[1]])

outb <- data.frame(seqnames(uniq.sites),start(uniq.sites),bnsm)

bnsm <- NULL
rm(bnsm)
gc()

colnames(out) <- c("chr","pos",sample.names)
rownames(out) <- NULL

colnames(outb) <- c("chr","pos",rep(sample.names,bootstrap.iter))
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

