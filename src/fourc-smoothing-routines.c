#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <stdio.h>

/**
 * idxes -- list of indexes for the matrix dm
 * dm -- NxM matrix of values
 * dim -- 2 integer vector of dm dimensions
 */
SEXP fourc_smoothing_mean(SEXP idxes, SEXP v) {
  int N = length(idxes);
  //int n = INTEGER(dim)[0];
  //int m = INTEGER(dim)[1];
  SEXP r = PROTECT(allocVector(REALSXP,N));
  //double *dmp = REAL(dm);
  double *pv = REAL(v);
  double *pr = REAL(r);
  
  for( int i = 0; i < N; i++ ) {
    SEXP vi = VECTOR_ELT(idxes,i);
    int *idx = INTEGER(vi);
    int l = length(vi);
    double d = 0;
    //double *cur = dmp;
    
    //for( int j = 0; j < m; j++ ) {
    //  double d = 0;
    //  if( j != 0 ) cur += n; 
      
    for( int k = 0; k < l; k++ ) {
      d += pv[idx[k]-1];
    }
      
    d /= l;
    pr[i] = d;
    //}
  }
  UNPROTECT(1);
  return r;
}

SEXP fast_write(SEXP fnamer,SEXP chrs, SEXP poses, SEXP m, SEXP dim) {
  FILE *fp;
  const char *fname = CHAR(STRING_ELT(fnamer,0));
  int *ppos = INTEGER(poses);
  double *pm = REAL(m);
  int l = length(chrs);
  int N = INTEGER(dim)[0];
  int M = INTEGER(dim)[1];
  
  if( ( fp = fopen(fname,"w")) == NULL) {
    error("Failed to open file for writing");
  }
  
  for( int i = 0; i < l; i++ ){
    const char *chr = CHAR(STRING_ELT(chrs,i));
    
    fprintf(fp,"%s\t%d\t",chr,ppos[i]);
    
    double *p = pm;
    
    for( int j = 0; j < M; j++ ) {
      if( j != 0 ) p += N;
      double d = p[i];
      if(j < (M-1)) {
        fprintf(fp,"%f\t",d); 
      } else {
        fprintf(fp,"%f\n",d);
      }
    }
  }
  
  
  fclose(fp);
  
  return R_NilValue;
}
