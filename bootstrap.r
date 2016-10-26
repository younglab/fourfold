library(compiler)
enableJIT(3)

NUM.OF.INTERATIONS <- 10

args <- commandArgs(T)

if(length(args)<2) {
  stop("Not enough arguments: need infile and outfile")
}

infile <- args[1]
outfile <- args[2]

counts <- read.table(infile,sep='\t')

N <- sum(counts$V3)
R <- nrow(counts)

#m <- matrix(0,nrow=R,ncol=NUM.OF.INTERATIONS)
prob <- counts$V3/N

#for( i in 1:NUM.OF.INTERATIONS) {
#  idx <- sample.int(R,size=N,replace=T,prob=prob)
  #n <- table(idx)
#  n <- tabulate(idx,nbins=R)
  
#  m[,i] <- as.vector(n) 
#  print(i)
#}

m <- do.call(cbind,lapply(1:NUM.OF.INTERATIONS,function(i) {
  print(i)
  
  idx <- sample.int(R,size=N,replace=T,prob=prob)
  tabulate(idx,nbins=R)
}))

df <- data.frame(counts[,1],counts[,2],m)

write.table(df,file=outfile,sep='\t',row.names=F,col.names=F,quote=F)


