#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

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
  
  for( int i = 0; i < N; i++ ) {
    SEXP idx = VECTOR_ELT(idxes,i);
    int l = length(idx);
    
    for( int j = 0; j < m; j++ ) {
      double d = 0;
      
      for( int k = 0; k < l; k++ ) {
        d += REAL(dm)[j*n+INTEGER(idx)[k]];
      }
      
      d /= l;
      REAL(r)[i*N+j] = d;
    }
  }
  UNPROTECT(1);
  return r;
}
