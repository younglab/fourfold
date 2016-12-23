#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>

double mean_helper(double *d,size_t n ) {
  double v = 0;
  
  for( int i = 0; i < n; i++ ) v += d[i];
  
  v /= n;
  return v;
}

//assumes this is a sample and not the population
double variance_helper(double *d,size_t n,const double *mean ) {
  double v = 0;
  double m = 0;
  
  
  if(mean==NULL) {
    for( int i = 0; i < n; i++ ) m += d[i];
    m /= n;
  } else {
    m = *mean;
  }
  
  for( int i = 0; i < n; i++ ) v += pow(d[i]-m,2);
  
  v /= (n-1);
  return v;
}

double median_helper(double *d,size_t n ) {
  double v = 0;
  
  int i = n/2;
  
  if( n % 2 == 0 ) { // even
    double x = d[i];
    double y = d[i+1];
    
    v = (x+y)/2.0;
  } else {
    v = d[i+1];
  }
  return v;
}

void quantile_helper(double *d,size_t n,double perclo, double perchi, double *l, double *h ) {
  int b = 0;
  
  //debugging
  //Rprintf(">");
  //for( int j = 0; j < n; j++ ) Rprintf("%f ",d[j]);
  //Rprintf("\n");

  for( int i = 0; i < n; i++ ) {
    double p = ((double)i)/(n-1);
    
    if(!b && p>=perclo) {
      *l = d[i];
      b = 1;
    }
    
    if(p>=perchi) {
      *h = d[i];
      break;
    }
  }
}

int comp(const void *a, const void *b) {
  double x = *((const double *)a);
  double y = *((const double *)b);
  
  if(x<y)
    return -1;
  else if(x>y)
    return 1;
  return 0;
}

/**
 * idxes -- list of indexes for the matrix dm
 * dm -- NxM matrix of values
 * dim -- 2 integer vector of dm dimensions
 */
SEXP fourc_smoothing_mean(SEXP fnamer,SEXP chrs, SEXP poses, SEXP idxes, SEXP dm, SEXP dim, SEXP statsv) {
  const char *fname = CHAR(STRING_ELT(fnamer,0));
  FILE *fp;
  int N = length(idxes);
  int n = INTEGER(dim)[0];
  int m = INTEGER(dim)[1];
  int *pposes = INTEGER(poses);
  int statsonly = INTEGER(statsv)[0]==1;

  double *pm = REAL(dm);
  
  if( ( fp = fopen(fname,"w") ) == NULL ) {
    error("Failed to open file");
  }
  
   //Rprintf("blah\n");
  
  for( int i = 0; i < N; i++ ) {
    SEXP vi = VECTOR_ELT(idxes,i);
    int *idx = INTEGER(vi);
    int l = length(vi);
    double d = 0;
    double *c = pm;
    double v[m];
    
    memset(&v,0,m);
    
    fprintf(fp,"%s\t%d\t",CHAR(STRING_ELT(chrs,i)),pposes[i]);
    
    for( int j = 0; j < m; j++ ) {
      double d = 0;
      if( j != 0 ) c += n; 
      
      for( int k = 0; k < l; k++ ) {
        d += c[idx[k]-1];
      }
      
      d /= l;
      
      if( statsonly ) {
        v[j] = d;
      } else {
        if( j == (m-1) ) 
          fprintf(fp,"%f\n",d);
        else
          fprintf(fp,"%f\t",d);
      }
    }
    
    if(statsonly) {
      //Rprintf("hmm\n");
      qsort(v,m,sizeof(v[0]),&comp);
      
      double me = mean_helper(v,m);
      double ve = variance_helper(v,m,&me);
      double md = median_helper(v,m);
      double lp = 0, hp = 0;
      quantile_helper(v,m,.025,.975,&lp,&hp);
      
      fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",me,ve,md,lp,hp);
    }
  }
  
  fclose(fp);
  
  return R_NilValue;
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
