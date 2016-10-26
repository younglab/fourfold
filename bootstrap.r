library(compiler)
enableJIT(3)

NUM.OF.INTERATIONS <- 1000

args <- commandArgs(T)

if(length(args)<2) {
  stop("Not enough arguments: need infile and outfile")
}

infile <- args[1]
outfile <- args[2]

counts <- read.table(infile,sep='\t')

N <- sum(counts$V3)
R <- nrow(counts)

prob <- counts$V3/N

m <- do.call(cbind,lapply(1:NUM.OF.INTERATIONS,function(i) {
  
  idx <- sample.int(R,size=N,replace=T,prob=prob)
  tabulate(idx,nbins=R)
}))

df <- data.frame(counts[,1],counts[,2],m)

write.table(df,file=outfile,sep='\t',row.names=F,col.names=F,quote=F)


