//Estimate expected utility for sequential boundaries
void euC(double *ustop, double *fdrstop, double *fnrstop, double *powerstop, double *jstop,double *b, int *ave,int *nrow,int *simid, double *j, double *u, double *fdr, double *fnr, double *power, double *summary, int *f,double *ineq, int *J, int *optlast);
void euCgrid(double *b0opt, double *b1opt, double *uopt, double *fdropt, double *fnropt, double *poweropt, double *jopt, double *ustop, double *fdrstop, double *fnrstop, double *powerstop, double *jstop, double *cf_sam, int *nb0, int *nb1, double *b0, double *b1, int *ave, int *nsim, int *nrow,int *simid, double *j, double *minj, double *u, double *fdr, double *fnr, double *power, double *summary, double *fdrmax, double *fnrmax, double *powmin, int *f, double *ineq, int *J, int *optlast);
void euCsearch(double *bopt, double *uopt, double *fdropt, double *fnropt, double *poweropt, double *jopt, double *cf_sam, int *npar, int *ngrid, double *binc, double *bini, int *nsim, int *nrow,int *simid, double *j, double *minj, double *u, double *fdr, double *fnr, double *power, double *summary, int *f, double *ineq, int *J, int *optlast, int *search, int *maxit);

//Forward Simulation
void forwsim_geneC(int *simid, int *time, double *u, double *fdr, double *fnr, double *power, double *summary, int *B, int *Bsummary, int *stepsum, int *J, int *Jini, int *mpred, int *util, double *cf, double *cf_sam, int *genelimit, double *v0thre, int *nrow, int *ncol, double *x, int *groups, int *K, double *alpha0, double *nu, double *balpha, double *nualpha,  int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *fdrmax, int *gapprox);
void forwsim_sampleC(int *simid, int *time, double *u, double *cc, double *summary, int *B, int *Bsummary, int *stepsum, int *J, int *Jini, int *mpred, int *genelimit, double *v0thre, int *nrow, int *ncol, double *x, int *groups, int *K, double *Kprob, double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, int *gapprox);

//Generalized Gamma-Gamma model fit
void fit_ggC(double *alpha0, double *nu, double *balpha, double *nualpha,  double *rho, double *prob, double *lhood, int *B, double *a_alpha0, double *b_alpha0, double *a_nu, double *b_nu, double *a_balpha, double *b_balpha, double *a_nualpha, double *b_nualpha, double *p_rho, double *p_prob, double *alpha0ini, double *nuini, double *balphaini, double *nualphaini, double *rhoini, double *probini, int *nrow, int *ncol, double *x, int *groups, int *K, int *nclust, int *npat, int *patterns, int *ngrouppat, int *gapprox, int *trace);

//Terminal utility routines
void maxwtpfp(double *u, int *d, double *fdr, double *fnr, double *cf, int *nsel, int *sel, double *v, int *npat);
void minfnrstfdr(double *u, double *threshold, int *d, double *fdr, double *fnr, double *power, int *nsel, int *sel, double *v, int *npat, double *fdrmax);
void sampleclas_ggC(int *d, double *pgroup, double *xnew, int *nsel, int *sel, int *nrow, int *ncol, double *x, int *groups, int *K, double *Kprob, double *rho, double *prob, double *alpha0, double *nu, double *balpha, double *nualpha,  int *nclust, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox);
void sel_mostDEgenes(int *nsel, int *sel, int *genelimit, double *v0thre, int *nrow, int *npat, double *v);
void uobsgeneC(double *u, double *fdr, double *fnr, int *util, int *nsel, int *sel, int *d, int *trued, double *cf);
void utgene_parC(double *u, int *d, double *fdr, double *fnr, double *power, double *threshold, int *util, double *cf, int *nsel, int *sel, double *v, int *npat, double *fdrmax);
void expected_fp(double *efp, double *threshold, int *nthre, int *B, int *niter, double *z, double *m, double *s, int *index, int *znclust, int *zclustsize, int *nrow, int *ncol, double *x, int *groups, int *K, double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, int *cluslist, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, double *sumxpred, double *prodxpred, double *nobsxpred, int *gapprox);
void bootnull(double *xboot, int *nrow, int *ncol, double *z, double *m, double *s, int *index, int *znclust, int *zclustsize);
void countde_threshold(int *nde, double *threshold, int *nthre, int *nrow, double *v, int *npat);
void countde(int *nde, double *threshold, int *nthre, double *fdrseq, int *nrow, double *v, int *npat);
void utsample_ggC(double *ccall, double *seccall, double *ccgroup, int *ngroup, int *B, double *preceps, int *genelimit, double *v0thre, int *nsel, int *sel, int *usesel, int *nrow, int *ncol, double *x, int *groups, double *v, int *K, double *Kprob, double *alpha0, double *nu, double *balpha, double *nualpha,  int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, int *ncolsumx, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox);

//Predictive terminal utility routines
void utgene_predC(double *m, double *s, int *deltat, int *B, double *preceps, int *util, double *cf, int *genelimit, double *v0thre, int *nsel, int *sel, int *usesel, int *nrow, int *ncol, double *x, int *groups, int *K, double *v, double *alpha0, double *nu, double *balpha, double *nualpha,  int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *fdrmax, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox);
void utsample_predC(double *ccall, double *seccall, double *ccgroup, int *ngroup, int *deltat, int *B, double *preceps, int *genelimit, double *v0thre, int *nsel, int *sel, int *usesel, int *nrow, int *ncol, double *x, int *groups, int *K, double *Kprob, double *v, double *alpha0, double *nu, double *balpha, double *nualpha,  int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox);

// Posterior computations
void compute_sumxC(double *sumx, double *prodx, double *nobsx, int *nsel, int *sel, int *ncolsumx, int *ncol, double *x, int *groups, int *K, int *npat, int *patterns, int *ngrouppat, int *init0);
void add_sumxC(double *sumxnew, double *prodxnew, double *nobsnew, int *nsel, int *sel, int *ncolsumx, double *sumx, double *prodx, double *nobsx);
void copy_sumxC(double *sumxnew, double *prodxnew, double *nobsnew, int *nsel, int *sel, int *ncolsumx, double *sumx, double *prodx, double *nobsx);
void pp_ggC(double *v, double *lhood, int *nsel, int *sel, int *ncol, double *x, int *groups, int *K, double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, int *cluslist, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, double *sumxpred, double *prodxpred, double *nobsxpred, int *usesumx, int *gapprox);
double pdfcond_pat_clus(int geneid, int patid, int clusid, double *alpha0, double *nu, double *balpha, double *nualpha, int *ngrouppat, int *colini, int ncolsumx, double *sumx, double *sumxpred, double *prodx, double *prodxpred, double *nobsx, double *nobsxpred, int usexpred);

// Simulate from the predictive
void simhyperpar_ggC(double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, double *rho, double *prob, int *cluslist, double *a_alpha0, double *b_alpha0, double *a_nu, double *b_nu, double *a_balpha, double *b_balpha, double *a_nualpha, double *b_nualpha, double *p_rho, double *p_prob, int *nrow, double *sumd, double *sumci, double *ngroupstot, double *sumalpha, double *sumlogalpha, double *suminvlambda, double *sumlambda, double *sumloglambda, int *npat, int *ngrouppat, int *gapprox);
void simnewsamples_ggC(double *xnew, int *dnew, double *anew, double *lnew, int *nsamples, int *groups, int *K, int *nsel, int *sel, double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, double *rho, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *gapprox);
void simpar_ggC(double *ngroupstot, double *sumd, double *sumci, double *sumalpha, double *sumlogalpha, double *suminvlambda, double *sumlambda, double *sumloglambda, int *groups, int *K, int *nrow, double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, double *rho, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *gapprox);
void simpred_ggC(double *xpred, int *d, double *alpha, double *l, int *usel, int *deltat, int *groups, int *K, int *nsel, int *sel, int *nrow, int *ncol, double *x, double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, double *rho, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox);
void simpred_oldggC(double *xpred, int *d, double *l, int *usel, int *deltat, int *groups, int *K, int *nrow, int *ncol, double *x, double *alpha, double *alpha0, double *nu, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *nobsx, int *usesumx);

// Gamma shape parameter conjugate prior
void dcgammaC(double *y, double *normk, double *x, int *n, double *a, double *b, double *c, double *d, double *r, double *s,  int *newton);
double kcgammaC(double *a, double *b, double *c, double *d, double *r, double *s, int *newton, int *log) ;
void mcgammaC(double *normk, double *m, double *v, double *a, double *b, double *c, double *d, double *r, double *s,  int *newton);
void rcgammaC(double *x, int *n, double *a, double *b, double *c, double *d, double *r, double *s,  int *newton);
void gapprox_par(double *aest, double *best, double *a, double *b, double *c, double *d, double *r, double *s, int *newton);
double logcgammaf(double x, double a, double b, double c, double d, double r, double s);
double logcgammafp(double x, double a, double b, double c, double d, double r, double s);
double logcgammafpp(double x, double a, double b, double c, double d, double r, double s);
