#include <R.h>                  /* for NA_REAL and Rmath.h */

/* This is the tolerance level */
const double EPSILON=0.0000001;

/* This function computes two-sample t-stat (without Welch correction) */
double two_sample_tstat(double *Y, const int *CL, const int *Ylen)
{
  double mean_no_na[2]={0,0},ss_no_na[2]={0,0},devi;
  double c0;
  int i,count[2]={0,0},class;
  double num, denom;   /* numerator and denominator */

  /*compute the mean and count first*/
  /* count is the number of objects in each class*/
  for (i=0; i<*Ylen; i++) {
    if (Y[i]==NA_REAL)          /* simply ignore NAs */
      continue;
    class=CL[i];
    mean_no_na[class] += Y[i];
    count[class]++;
  }
  mean_no_na[0]/=(count[0]*1.0);
  mean_no_na[1]/=(count[1]*1.0);

  /*compute the variance in each group*/
  for (i=0; i<*Ylen; i++) {
    if (Y[i]==NA_REAL)
      continue;
    class=CL[i];
    devi=(Y[i]-mean_no_na[class]);
    ss_no_na[class] += devi*devi;
  }

  /* This is the equal-var t-stat; not the Welch-corrected stat. Note that the sign is compatible with genefilter/t.test(), not multtest. */
  /* num=mean_no_na[1]-mean_no_na[0]; */
  num=mean_no_na[0]-mean_no_na[1];
  c0= (double)(count[0]+count[1])/(count[0]*count[1]*(count[0]+count[1]-2));
  denom=sqrt((ss_no_na[0]+ss_no_na[1])*c0);
  /* check the denom. */
  if(denom<EPSILON)
    return NA_REAL; /*If too small denominator, treat as NA*/
  return num/denom;
}

/* This function computes two-sample t-stat (with Welch correction) */
double two_sample_tstatWel(double *Y, const int *CL, const int *Ylen)
{
  double mean_no_na[2]={0,0},ss_no_na[2]={0,0},devi;
  /*double df;*/
  int i,count[2]={0,0},class;
  double num, denom;   /* numerator and denominator */

  /*compute the mean and count first*/
  /* count is the number of objects in each class*/
  for (i=0; i<*Ylen; i++) {
    if (Y[i]==NA_REAL)          /* simply ignore NAs */
      continue;
    class=CL[i];
    mean_no_na[class] += Y[i];
    count[class]++;
  }
  mean_no_na[0]/=(count[0]*1.0);
  mean_no_na[1]/=(count[1]*1.0);

  /*compute the variance in each group*/
  for (i=0; i<*Ylen; i++) {
    if (Y[i]==NA_REAL)
      continue;
    class=CL[i];
    devi=(Y[i]-mean_no_na[class]);
    ss_no_na[class] += devi*devi;
  }

  /* This is the Welch-corrected t-stat. Note that the sign is compatible with genefilter/t.test(), not multtest. */
  /* num=mean_no_na[1]-mean_no_na[0]; */
  num=mean_no_na[0]-mean_no_na[1];
  /* c0= (double)(count[0]+count[1])/(count[0]*count[1]*(count[0]+count[1]-2)); */
  denom=sqrt(ss_no_na[0]/(count[0]-1.0)/count[0]+ss_no_na[1]/(count[1]-1.0)/count[1]);
  /* check the denom. */
  if(denom<EPSILON)
    return NA_REAL; /*If too small denominator, treat as NA*/
  return num/denom;
}

/* compute delta (difference) from two genes Y1, Y2 and then
   return two-sample t-statistic. */
double delta_tstat(const double *Y1, const double *Y2, const int *nslides, const int *CL)
{
  /* assign a vector of size *nslides to hold delta */
  double delta[*nslides];
  double t;
  int i;
  for(i=0; i<*nslides; i++){
    delta[i]=Y1[i]-Y2[i];
  }
  /* Now compute the t-stats */
  t=two_sample_tstatWel(delta, CL, nslides);
  return t;
}

/* Function superdelta_stats_ttest() below takes gene expression data (Y), a vector of indices of baseline genes (basegenes),
 a vector of group class labels (CL), and a pre-defined matrix Tmat. Note that Y must be organized in such a way that:
 1. The first nbase genes are the baseline genes;
 2. CL must be a 0,1 vector;
 3. Tmat must be a ngenes * nbase matrix. */
void superdelta_tstats(const double *Y, const int *ngenes, const int *nslides, const int *nbase, const int *CL, double *Tmat)
{
  int n=*nslides, p=*ngenes;
  int i,j;
  /* First, let's copy the vectorized Y into 2d-array in C (Yc). */
  double **Yc;
  Yc = (double **) R_alloc(p, sizeof(double));
  for(i = 0; i < p; i++) {
    Yc[i] = (double *) R_alloc(n, sizeof(double));
    for(j = 0; j < n; j++){
      Yc[i][j] = Y[j*p+i];
    }
  }

  /* Now let's compute the (nbase-1) x (nbase-1) dim submatrix of
     Tmat, which is symmetric. */
  for(i=0; i<*nbase; i++){
    /* I decide to assign NA to diagonal elements. */
    Tmat[i*p+i]=NA_REAL;
    /* Off-diagonal elements; use symmetry to save time */
    for(j=0;j<i;j++){
      Tmat[j*p+i]=delta_tstat(Yc[i], Yc[j], nslides, CL);
      Tmat[i*p+j]= -Tmat[j*p+i];
    }
  }
  /* The rest of the Tmat; no symmetry. */
  for(i=*nbase; i<p; i++){
    for(j=0;j<*nbase;j++){
      Tmat[j*p+i]= delta_tstat(Yc[i], Yc[j], nslides, CL);
    }
  }
}

/* A convenient function that can be used in place of rowttest() in genefilter package */
void fastt(const double *Y, const int *ngenes, const int *nslides, const int *CL, double *Tvec)
{
  int n=*nslides, p=*ngenes;
  int i,j;
  /* First, let's copy the vectorized Y into 2d-array in C (Yc). */
  double **Yc;
  Yc = (double **) R_alloc(p, sizeof(double));
  for(i = 0; i < p; i++) {
    Yc[i] = (double *) R_alloc(n, sizeof(double));
    for(j = 0; j < n; j++){
      Yc[i][j] = Y[j*p+i];
    }
  }

  /* Now let's compute the two-sample t-statistics (with Welch adjustment). */
  for(i=0; i<*ngenes; i++){
    Tvec[i]=two_sample_tstatWel(Yc[i], CL, nslides);
  }
}
