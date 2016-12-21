#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>

/**
 * idxes -- list of indexes for the matrix dm
 * dm -- NxM matrix of values
 * dim -- 2 integer vector of dm dimensions
 */
SEXP fourc_smoothing_mean(SEXP idxes, SEXP dm, SEXP dim) {
  int N = length(idxes);
  int n = INTEGER(dim)[0];
  int m = INTEGER(dim)[1];
  SEXP r = PROTECT(allocMatrix(REALSXP,N,m));
  double *dmp = REAL(dm);
  
  for( int i = 0; i < N; i++ ) {
    SEXP idx = VECTOR_ELT(idxes,i);
    int l = length(idx);
    double *cur = dmp;
    
    for( int j = 0; j < m; j++ ) {
      double d = 0;
      if( j != 0 ) cur += n; 
      
      for( int k = 0; k < l; k++ ) {
        //int t = j*n+(INTEGER(idx)[k]-1);
        
        //Rprintf("[%d,%d,%d,%d,%d,%d,%d,%d]\n",i,j,k,INTEGER(idx)[k],N,n,m,t);
        
        //d += dmp[j*n+(INTEGER(idx)[k]-1)];
        d += cur[INTEGER(idx)[k]-1];
      }
      
      d /= l;
      //Rprintf("----[%d,%d] %d %d %d\n",i,j,N,n,m);
      REAL(r)[i+j*N] = d;
    }
  }
  UNPROTECT(1);
  return r;
}
