/***********************************************************
 Basic statistical, input-output and matrix manipulation

 Authors. Peter Mueller, Stephen Morris, David Rossell
 Last modified. 04 2007
***********************************************************/

#if !defined(STAT_H_INCLUDED)
#define STAT_H_INCLUDED

#include <math.h>
#include <stdio.h>

#if !defined(M_PI)
#define M_PI (3.1415926535897932385)
#endif

#if !defined(M_PI_2)
#define M_PI_2 (1.570796326794896619231321691640)
#endif

#if !defined(SQ_M_PI_2)
#define SQ_M_PI_2 (2.5066282746310002416123552)
#endif

/**************************************************************/
/* Functions to compute means & variances                     */
/**************************************************************/

double meani(int *x, int lim);
double wmeani(int *x, int lim, double *w);
double meanx(double *x, int lim);
double wmeanx(double *x, int lim, double *w);
void colMeans(double *m, double *x, int nrow, int ncol);
double vari(int *x, int lim, int unb);
double wvari(int *x, int lim, double *w);
double varx(double *x, int lim, int unb);
double wvarx(double *x, int lim, double *w);

/**************************************************************/
/* Input/output functions (interface)                         */
/**************************************************************/

FILE *openIn(char *);
FILE *openOut(char *);
void scanFloat(char *, float *);
void scanDouble(char *, double *);
void scanInt(char *, int *);
void fscanDouble(FILE *,char *, double *);
void fscanInt(FILE *,char *, int *);
void scanLong(char *, long *);

void scanFloatArray(char *,float *, int);
void scanArray(char *,float *, int);
void scanDoubleArray(char *,double *, int);
void scanString(char *txt, char *s, int n);
void fscanString(FILE *,char *txt, char *s, int n);
void fscanDoubleArray(FILE *,double *, int);
void scanDoubleMatrix(char *, double **, int, int);
void fscanDoubleMatrix(FILE *ifile, double **x,int r,int c);
void scanIntArray(char *,int *, int);
void fscanIntArray(FILE *ifile, int *x, int n);

void writeInt(int);
void writeLong(long i);
void writeFloat(float);
void writeDouble(double);

void writeIntArray(int *,int, int);
void fwriteIntArray(FILE *, int *,int, int);
void fwriteIntMatrix(FILE *f, int **x, int rows, int cols);
void writeIntMatrix(int **x, int rows, int cols);
void writeDoubleArray(double *,int, int);
void writeDoubleMatrix2(double **, int , int);
void fwriteDoubleArray(FILE *, double *,int, int);
void fwriteDoubleMatrix2(FILE *, double **, int , int);
void writeDoubleMatrix(double **, int, int);
void writeFloatArray(float *, int, int);
void writeArray(float *, int, int); 

void fserror(char *proc, char *act, char *what);

/**************************************************************/
/* Debug messages etc. (mess)                                 */
/**************************************************************/

void error(char *,char *, int);
void err_msg(char *fct, char *txt, int n1, int n2, int n3);

/**************************************************************/
/* Memory allocation                                          */
/**************************************************************/

float   *vector(int,int);
double  *dvector(int,int);
double  **dmatrix(int,int,int,int);
double  ***darray_3(int, int);
double ***darray3(int n, int p, int q);
int     *ivector(int,int);
int     **imatrix(int,int,int,int);
int ***iarray_3(int lo, int hi);
int ***iarray3(int p1, int p2, int p3);

void free_vector(float  *,int,int);
void free_dvector(double  *,int,int);
void free_ivector(int  *,int,int);
void free_dmatrix(double  **,int,int,int,int);
void free_imatrix(int  **,int,int,int,int);

void nrerror(char *proc, char *act, char *what);

/**************************************************************/
/* Mathematical functions                                     */
/**************************************************************/

double digamma(double x);                          //from S Poetry (by Patrick J. Burns)
double trigamma(double x);
double polygamma(double x, long n, double low, double high, long terms, double nfact); //from S Poetry
double lnbeta(double a, double b); //log of Beta function
double betacf(double a, double b, double x); //continued fraction for incomplete Beta function

/**************************************************************/
/* Vector algebra (vector)                                    */
/**************************************************************/

void grid (double x0, double xn, int n, double *x);
void rA(double,double **, double **, int, int);  //Multiply matrix by scalar
void Ax_plus_y(double **A, double *x, double *y, double *z, int ini, int fi); //matrix*vector+vector
void xA(double *x,double **A,double *z, int ini, int fi);  //Multiply vector * matrix
void Ax(double **A,double *x,double *z, int ini, int fi);  //Multiply matrix * vector
void a_plus_b(double *a, double *b, double *c, int ini, int fi); //Vector sum i.e. c[i]=a[i]+b[i]
void a_prod_b(double *a, double *b, double *c, int ini, int fi); //Vector prod i.e. c[i]=a[i]*b[i]
void R_zero(double **, int, int);
int iabs(int x);   //absolute value of an integer
int imax_xy(int x, int y);
int imin_xy(int x, int y);
double max_xy(double x, double y);
double min_xy(double x, double y);

void choldc(double **a, int n, double **aout);   //Cholesky decomposition
void choldc_inv(double **a, int n, double **aout); //Inverse of Cholesky decomposition
double choldc_det(double **chols, int n); //Determinant of a symmetric def+ using its Cholesky decomp
void inv_posdef(double **a, int n, double **aout); //Inverse of a positive definite matrix

void ludc(double **a, int n, int *indx, double *d); //LU decomposition (renamed routine ludcmp from NR)
void lu_solve(double **a, int n, int *indx, double b[]); //Solve A*x=b (renamed routine lubksb from NR)
void lu_inverse(double **a, int n, double **aout); //Inverse of A[1..n][1..n]
double lu_det(double **a, int n); //Determinant of A[1..n][1..n]

int dcompare(const void *a, const void *b);               
void dvecsort(double *v, int size);                           //sort a vector using qsort from stdlib
void dindexsort(double *x, int *index, int ilo, int ihi, int incr); //sort a vector of indexes using self-written quicksort routine
void iindexsort(int *x, int *index, int ilo, int ihi, int incr); //like dindexsort but for integers


/**************************************************************/
/* Random sampling                                            */
/**************************************************************/

void samplei_wr(int *x, int popsize, int n); //sample wo replacement from a vector of integers
void sampled_wr(double *x, int popsize, int n); //same for vector of doubles

/**************************************************************/
/* Random variate generation (rand)                           */
/**************************************************************/

//Several
void setseed(long, long);
double  runif();
int runifdisc(int min, int max);
int rdisc(double *probs, int nvals);
double rbeta(double , double );
double pbeta(double x, double pin, double qin); //quantile from a Beta(pin,qin)
void rdirichlet(double *w, double *alpha, int *p);
double gamdev(double );
int	binomial(int , double );
void multinomial(int, int, double *, int *);

// Normal
double dnorm(double y, double m, double s, int logscale); //density of Normal(m,s^2)
double dmvnorm(double *y, int n, double *mu, double **cholsinv, double det, int logscale); //density of multivariate Normal
double	qnorm (double cdf, double m, double s);  //quantile from Normal(m,s^2)
double	pnorm(double y, double m, double s);  //cdf of Normal(m,s^2)
double rnorm(double mu, double s); //draw from univariate Normal(mu,s^2)
double rnorm_trunc(double ltrunc, double rtrunc, double m, double s); //draw trunc Normal given trunc points
double rnorm_trunc_prob(double lprob, double rprob, double m, double s); //draw trunc Normal given trunc probs
void rmvnorm(double *y, int n, double *mu, double **chols); //draw from multivariate Normal

// T Student
double dt(double y, double mu, double s, int nu); //density of t with nu df
double dmvt(double *y, int n, double *mu, double **cholsinv, double det, int nu, int logscale); //density of multivariate t
double rt(int nu); //draw from univariate t with nu degrees of freedom
double rt_trunc(int nu, double ltrunc, double rtrunc); //draw from truncated t given trunc points
double rt_trunc_prob(int nu, double lprob, double rprob);  //draw from truncated t given trunc probs
double qt(double p, int nu);  //quantile from t-Student with nu degrees of freedom
double pt(double x, int nu);  //CDF of t-Student with nu degrees of freedom
void rmvt(double *y, int n, double *mu, double **chols, int nu); //draw from multivar T with nu degrees of freedom

/* More random variate stuff (dcdflib, from CMU statlib "www.stat.cmu.edu") */
double fifdint(double);
void cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
double spmpar(int*);
void cumnor(double*,double*,double*);
double dinvnr(double *p,double *q);
double stvaln(double*);
double devlpl(double [],int*,double*);
double gamln(double*);
double gamln1(double*);
extern int ipmpar(int*);                      /* code in ipmpar.c */

/*even more stuff (ranlib) */
extern double genunf(double low,double high);
extern double gengam(double a,double r);
extern double sgamma(double a);
extern double snorm(void);
double fsign( double num, double sign );
extern double sexpo(void);
extern long mltmod(long a,long s,long m);
extern double ranf(void);
extern void gscgn(long getset,long *g);
extern void setall(long iseed1,long iseed2);   /* code in com.c */
extern void initgn(long isdtyp);               /* code in com.c */
extern long ignlgi(void);                      /* code in com.c */
extern void inrgcm(void);                      /* code in com.c */


/**************************************************************/
/* Integration                                                */
/**************************************************************/

double midpnt(double (*func)(double), double a, double b, int n); //nth stage refinement of integral of func from a to b (evenly spaced in x)
double midinf(double (*funk)(double), double aa, double bb, int n); //nth stage refinement of integral of funk from aa to bb (evenly spaced in 1/x)
double qromo(double (*func)(double), double a, double b, double (*choose)(double(*)(double), double, double, int)); //Romberg integr on open interval (a,b)

/**************************************************************/
/* Interpolation and extrapolation                            */
/**************************************************************/

void polint (double xa[], double ya[], int n, double x, double *y, double *dy); //interpolates via polynomials


/**************************************************************/
/* Function optimization                                      */
/**************************************************************/

double brent(double ax,double bx,double cx,double (*f)(double),double tol,double *xmin,int itmax); //univariate minim
double dbrent(double ax,double bx,double cx,double (*f)(double),double (*df)(double),double tol,double *xmin,int itmax);
void powell(double p[],double **xi,int n,double ftol,int *iter,double *fret,double (*func)(double []),int itmax);//multivar minim
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []), int itmax); //minim in 1 direction
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double (*func)(double)); //find bracketing triplets

#endif
