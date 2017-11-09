suppressMessages(library(compiler))
suppressMessages(library(Rcpp))
suppressMessages(library(bigmemory))


cppFunction('IntegerMatrix numericBootstrap(NumericVector probs, IntegerVector nSites, IntegerVector nReads, IntegerVector nBoot) {
  IntegerMatrix counts(nSites[0],nBoot[0]);
  int nsites = nSites[0];
  int nreads = nReads[0];
  int nboot = nBoot[0];

  for( int i = 0; i < nboot; i++ ) {
    std::vector<double> pv;
    pv.reserve(nreads);
    for( int x = 0; x < nreads; x++ ) pv[x] = R::runif(0,1);
    std::sort(pv.begin(),pv.end());
    int hit = 0;


    for( double p : pv ) {
      while( hit < nsites && probs[hit] <= p) hit++;

      counts(hit,i)++;
    }
  }

  return counts;
}',
plugins=c("cpp11"),
depends=c("BH","bigmemory"),
includes=c('#include "bigmemory/BigMatrix.h"','#include "bigmemory/MatrixAccessor.hpp"')
)


enableJIT(3)

args <- commandArgs(T)

if(length(args) < 5 ) {
  stop("Not enough arguments:<infile> <outfile counts> <outfile rpm> <number of mapped reads> <bootstrap iterations>")
}

infile <- args[1]
outfilecounts <- args[2]
outfilerpm <- args[3]
mappedreads <- as.integer(args[4])
bootstrapiter <- as.integer(args[5])

### Matrix structure
### Column 1: Chromosome name
### Column 2: Position of fragment/digestion
### Column 3: Read count/rpm

counts <- read.table(infile,sep='\t')

N <- sum(counts[,3])
R <- nrow(counts)

probs <- cumsum(counts[,3]/N)

m <- numericBootstrap(probs,R,N,bootstrapiter)

df <- data.frame(counts[,1],counts[,2],m)
dfrpm <- data.frame(counts[,1],counts[,2],m/mappedreads*1e6)

write.table(df,file=outfilecounts,sep='\t',row.names=F,col.names=F,quote=F)
write.table(dfrpm,file=outfilerpm,sep='\t',row.names=F,col.names=F,quote=F)


