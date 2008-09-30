#include <math.h>

#include <stdio.h>



#include "cstat.h"

#include "cseqdesma.h"





//int __iglobal=0;  //debug

double meanu, meanj;                      //Global variables used by euC, euCgrid

struct pars {

  double *probest;

  double *rhoest;

  double *lhood;

  int *justlhood; 

  int *nsel;

  int *sel;

  int *K;

  int *equalcv;

  int *nclust;

  double *rho;

  double *prob;

  int *npat;

  int *patterns;

  int *ngrouppat;

  int *colini;

  int *ncolsumx;

  double *sumx;

  double *prodx;

  double *nobsx;

  int *gapprox;

} esteppars;





/*********************************************************************************************

                      ESTIMATE EXPECTED UTILITY FOR SEQUENTIAL BOUNDARIES

*********************************************************************************************/



void euC(double *ustop, double *fdrstop, double *fnrstop, double *powerstop, double *jstop,double *b, int *ave,int *nrow,int *simid, double *j, double *u, double *fdr, double *fnr, double *power, double *summary, int *f,double *ineq, int *J, int *optlast)

{

/* Evaluates realized utility for a single parametric or non-parametric boundary */

/* Input arguments

  - b: boundary at which to evaluate the observed terminal utility, fnr, fdr, stopping time

  - ave: ave==0 means returning utilities and stopping times for each simulation. ave==1 means returning averages only.

  - nrow: length of simid, j, u, fdr, fnr, summary

  - simid: simulation identifier

  - j: time

  - u: expected terminal utility in simulation simid if we were to stop at time j

  - fdr: expected fdr in simulation simid if we were to stop at time j

  - fnr: expected fnr in simulation simid if we were to stop at time j

  - power: power in simulation simid if we were to stop at time j

  - summary: summary statistic used to take stopping decision

  - f: f==0 means non-parametric boundaries. f==1 means linear boundaries. f==2 means inverse square root boundaries. If f==0 b must be of length equal to the number of possible stopping times, if f!=0 b must be of length 2.

  - ineq: ineq==1 means stop if summary>boundary; ineq==-1 means stop if summary<boundary

  - J: number of distinct stopping time points

  - optlast: optlast==1 means continue at j=J-1 iff summary>0 (corresponds to optimal rule for several utilities)

*/

/* Output arguments

  - ustop: If ave==0 a vector with the observed terminal utility for each simulation. If ave==1 only the mean is returned.

  - fdrstop: fdr averaged over simulations.

  - fnrstop: fnr averaged over simulations.

  - powerstop: power averaged over simulations.

  - jstop: If ave==0 a vector with the stopping time for each simulation. If ave==1 only the mean is returned.

 */

								

  int i, k, jdif;

  int found= 0;



  k= jdif= 0; meanu= meanj= 0; (*fdrstop)= (*fnrstop)= (*powerstop)= 0;

  for(i=0;i<((*nrow)-1); i++) {

    if ((found==0) && ((jdif<(*J-2)) || ((jdif==(*J-2)) && ((*optlast)==0)) )) {     //If not in 2 last time points OR in 2nd last with *optlast==0

      if ((*f)==1) {                                                  //Linear boundaries

        if ((*ineq)*summary[i] > (*ineq)*(b[0]+b[1] * j[i])) {

          if ((*ave)==0) { ustop[k]= u[i]; jstop[k]= j[i]; }

          meanu += u[i]; meanj += j[i]; found= 1; (*fdrstop) += fdr[i]; (*fnrstop) += fnr[i]; (*powerstop) += power[i]; k++;   

        }

      } else if ((*f)==2) {                                           //Inv sqrt boundaries

        if ((*ineq)*summary[i] > (*ineq)*(b[0]+b[1]/sqrt(j[i]))) {

          if ((*ave)==0) { ustop[k]= u[i]; jstop[k]= j[i]; }

	  meanu += u[i]; meanj += j[i]; found= 1; (*fdrstop) += fdr[i]; (*fnrstop) += fnr[i]; (*powerstop) += power[i]; k++;  

        }

      } else if ((*f)==0) {                                           //Non-parametric boundaries

        if ((*ineq)*summary[i] > (*ineq)*b[jdif]) {

          if ((*ave)==0) { ustop[k]= u[i]; jstop[k]= j[i]; }

	  meanu += u[i]; meanj += j[i]; found= 1; (*fdrstop) += fdr[i]; (*fnrstop) += fnr[i]; (*powerstop) += power[i]; k++;  

        }

      }

    } else if ((found==0) && (jdif==(*J-2)) && ((*optlast)==1)) {                 //If at the second to the last time point and *optlast==1

      if (summary[i]<0) {

        if ((*ave)==0) { ustop[k]= u[i]; jstop[k]= j[i]; }

	meanu += u[i]; meanj += j[i]; found= 1; (*fdrstop) += fdr[i]; (*fnrstop) += fnr[i]; (*powerstop) += power[i]; k++;  

      }

    } else if ((found==0) && (jdif==(*J-1))) {

      if ((*ave)==0) { ustop[k]= u[i]; jstop[k]= j[i]; }                            //Always stop at the time horizon

      meanu += u[i]; meanj += j[i]; (*fdrstop) += fdr[i]; (*fnrstop) += fnr[i]; (*powerstop) += power[i]; k++;  

    }

    

    jdif++; 

    if (simid[i]!=simid[i+1]) { found = 0; jdif=0; }                                //Set found,jdif=0 for a new simulation 

  }   //End i for



  if (found==0) {

      if ((*ave)==0) { ustop[k]= u[i]; jstop[k]= j[i]; }

      meanu += u[i]; meanj += j[i]; (*fdrstop) += fdr[i]; (*fnrstop) += fnr[i]; (*powerstop) += power[i]; k++;  

  }



  meanu = meanu/k; meanj = meanj/k;                                                 //Obtain mean utility & stopping time

  if ((*ave)==1) { (*ustop)= meanu; (*jstop)= meanj; }    

  (*fdrstop) = (*fdrstop) / k; (*fnrstop) = (*fnrstop) / k; (*powerstop) = (*powerstop) / k;   //Obtain mean FDR, FNR & power

  

}





void euCgrid(double *b0opt, double *b1opt, double *uopt, double *fdropt, double *fnropt, double *poweropt, double *jopt, double *ustop, double *fdrstop, double *fnrstop, double *powerstop, double *jstop, double *cf_sam, int *nb0, int *nb1, double *b0, double *b1, int *ave, int *nsim, int *nrow,int *simid, double *j, double *minj, double *u, double *fdr, double *fnr, double *power, double *summary, double *fdrmax, double *fnrmax, double *powmin, int *f, double *ineq, int *J, int *optlast) {

/* Evaluates realized utility on a grid */

/* Input arguments

  - cf_sam: sampling cost

  - nb0: length of b0

  - nb1: length of b1

  - b0: vector of b0 values

  - b1: vector of b1 values

  - ave: ave==0 means returning utilities and stopping times for each simulation. ave==1 means returning averages only.

  - nsim: number of simulations (i.e. number of distinct elements in simid)

  - nrow: length of simid, j, u, fdr, fnr, summary

  - simid: simulation identifier

  - j: time

  - minj: minimum of j

  - u: terminal utility in simulation simid if we were to stop at time j

  - fdr: fdr in simulation simid if we were to stop at time j

  - fnr: fnr in simulation simid if we were to stop at time j

  - power: power in simulation simid if we were to stop at time j

  - summary: summary statistic used to take stopping decision

  - fdrmax: contraint on fdr. Optimization chooses b0,b1 satisfying fdr<=fdrmax (if such b0,b1 exists)

  - fnrmax: contraint on fnr. Optimization chooses b0,b1 satisfying fnr<=fnrmax (if such b0,b1 exists)

  - powmin: contraint on power. Optimization chooses b0,b1 satisfying power>=powermin (if such b0,b1 exists)

  - f: f==1 means linear boundaries. f==2 means inverse square root boundaries

  - ineq: ineq==1 means stop if summary>boundary; ineq==-1 means stop if summary<boundary

  - J: number of distinct stopping time points

  - optlast: optlast==1 means continue at j=J-1 iff summary>0 (corresponds to optimal rule for several utilities)

*/

/* Output arguments

  - b0opt: value of b0 maximizing the mean utility (naive estimate) amongst those satisfying fdr<=fdrmax,fnr<=fnrmax,power>=powermin

  - b1opt: value of b1 maximizing the mean utility (naive estimate) amongst those satisfying fdr<=fdrmax,fnr<=fnrmax,power>=powermin

  - uopt: mean utility at b0opt,b1opt

  - fdropt: mean fdr at b0opt,b1opt

  - fnropt: mean fnr at b0opt,b1opt

  - poweropt: mean power at b0opt,b1opt

  - jopt: mean j at b0opt,b1opt

  - ustop: if ave==0 a vector with the observed terminal utility for each simulation. If ave==1 only the means returned.

  - fdrstop: fdr averaged over simulations.

  - fnrstop: fnr averaged over simulations.

  - jstop: if ave==0 a vector with the stopping time for each simulation. If ave==1 only the means are returned.

 */

								

  int i, k, pos1, pos2, feasible=0, feasibleopt=0;

  double b[2];



  for(i=0;i<(*nb0); i++) {

      for (k=0; k<(*nb1); k++) {

        b[0] = b0[i]; b[1] = b1[k];

        pos2= i*(*nb1)+k;

        if ((*ave)==1) { pos1= pos2; } else { pos1= (*nsim)*pos2; }

        euC(ustop+pos1,fdrstop+pos2,fnrstop+pos2,powerstop+pos2,jstop+pos1,b,ave,nrow,simid,j,u,fdr,fnr,power,summary,f,ineq,J,optlast); //Evaluate utility

	feasible= (*(fdrstop+pos2) <= *fdrmax) + (*(fnrstop+pos2) <= *fnrmax) + (*(powerstop+pos2) >= *powmin);

        if ((i==0) && (k==0)) {                                                                               //Initialize optimal values

          (*b0opt) = b[0]; (*b1opt) = b[1];

          (*uopt)= meanu - (*cf_sam)*((*jopt)-(*minj)+1); (*jopt) = meanj;

          (*fdropt) = *(fdrstop+pos2); (*fnropt) = *(fnrstop+pos2); (*poweropt) = *(powerstop+pos2);

	  feasibleopt= feasible;

        }

        //If current point satisfies more constraints than the optimal, or it satisfies as many constraints AND it improves utility, update optimum 

        if ((feasible>feasibleopt) || ((feasible==feasibleopt) && (meanu-(*cf_sam)*(meanj-(*minj)+1)>(*uopt)))) {

          (*b0opt) = b[0]; (*b1opt) = b[1];

          (*uopt)= meanu - (*cf_sam)*((*jopt)-(*minj)+1); (*jopt) = meanj;

          (*fdropt) = *(fdrstop+pos2); (*fnropt) = *(fnrstop+pos2); (*poweropt) = *(powerstop+pos2);

	  feasibleopt= feasible;

        }        

      }                                                                               //End k for

  }                                                                                   //End i for



}





void euCsearch(double *bopt, double *uopt, double *fdropt, double *fnropt, double *poweropt, double *jopt, double *cf_sam, int *npar, int *ngrid, double *binc, double *bini, int *nsim, int *nrow,int *simid, double *j, double *minj, double *u, double *fdr, double *fnr, double *power, double *summary, int *f, double *ineq, int *J, int *optlast, int *search, int *maxit) {

/* Searches for the optimal parametric or non-parametric boundary by successive univariate optimizations */

/* Input arguments

  - cf_sam: sampling cost

  - npar: number of parameters i.e. length of bopt and bini

  - ngrid: in the univariate optim the exp utility is evaluated for a univariate grid of ngrid equidistant points (ngrid increased to nearest odd number)

  - binc: vector of length npar indicating how far apart must the grid points be from each other for each parameter

  - bini: vector of b values

  - nsim: number of simulations (i.e. number of distinct elements in simid)

  - nrow: length of simid, j, u, fdr, fnr, summary

  - simid: simulation identifier

  - j: time

  - minj: minimum of j

  - u: terminal utility in simulation simid if we were to stop at time j

  - fdr: fdr in simulation simid if we were to stop at time j

  - fnr: fnr in simulation simid if we were to stop at time j

  - power: power in simulation simid if we were to stop at time j

  - summary: summary statistic used to take stopping decision

  - f: f==0 means non-parametric boundaries. f==1 means linear boundaries. f==2 means inverse square root boundaries. If f==0 b must be of length equal to the number of possible stopping times, if f!=0 b must be of length 2.

  - ineq: ineq==1 means stop if summary>boundary; ineq==-1 means stop if summary<boundary

  - J: number of distinct stopping time points

  - optlast: optlast==1 means continue at j=J-1 iff summary>0 (corresponds to optimal rule for several utilities)

  - search: search==0 indicates to estimate exp utility at bini and not perform any optimization. search==1 indicates to perform optimization.

  - maxit: maximum number of iterations in the optimization scheme

*/

/* Output arguments

  - bopt: value of b maximizing the mean utility (naive estimate)

  - uopt: mean utility at b0opt,b1opt

  - fdropt: mean fdr at b0opt,b1opt

  - fnropt: mean fnr at b0opt,b1opt

  - poweropt: mean power at b0opt,b1opt

  - jopt: mean j at b0opt,b1opt

 */



  int one=1, i, k, found, it;

  double *bcur, *bmin, ucur, fdrcur, fnrcur, powercur, jcur;



  if ((*ngrid)%2 == 0) { (*ngrid) += 1; }                                                         //force *ngrid to be an odd number

  bcur= dvector(0,*npar); bmin= dvector(0,*npar);

  for (i=0; i<(*npar); i++) { bcur[i]= bini[i]; bopt[i]= bini[i]; }                               //initialize bopt, bcur

  euC(uopt,fdropt,fnropt,poweropt,jopt,bopt,&one,nrow,simid,j,u,fdr,fnr,power,summary,f,ineq,J,optlast);                   //initialize uopt



  if ((*search)==1) { found=0; it= 0; } else { found=1; }

  while ((found==0) && (it< *maxit)) {     //apply optimization scheme until local maximum found or reached limit of iterations

    found= 1;

    for (i=0; i<(*npar); i++) {          //for each parameter perform univariate optimization fixing the remaining parameters

      bmin[i]= bcur[i] - (*ngrid/2)*binc[i];

      for (k=0; k<(*ngrid); k++) {       //loop over 1-dim grid points

	bcur[i]= bmin[i] + k*binc[i];

	euC(&ucur,&fdrcur,&fnrcur,&powercur,&jcur,bcur,&one,nrow,simid,j,u,fdr,fnr,power,summary,f,ineq,J,optlast);

	if ((ucur-(*cf_sam)*jcur)>((*uopt)-(*cf_sam)*(*jopt))) { 

	  found=0; bopt[i]=bcur[i]; (*uopt)=ucur; (*fdropt)=fdrcur; (*fnropt)=fnrcur; (*poweropt)=powercur; (*jopt)=jcur;

	}

      }  //end loop over 1-dim grid points

      bcur[i]= bopt[i];                  //fix ith parameter to its 1-dim optimal value

    }  //end for each param perform univariate optimization

    it++;

  } //end apply optimization scheme



  (*uopt)= (*uopt) - (*cf_sam)*((*jopt)-(*minj)+1);    //substract sampling cost from utility



  free_dvector(bcur, 0,*npar); free_dvector(bmin, 0,*npar);



}







/*********************************************************************************************

                                    FORWARD SIMULATION

*********************************************************************************************/



void forwsim_geneC(int *simid, int *time, double *u, double *fdr, double *fnr, double *power, double *summary, int *B, int *Bsummary, int *stepsum, int *J, int *Jini, int *mpred, int *util, double *cf, double *cf_sam, int *genelimit, double *v0thre, int *nrow, int *ncol, double *x, int *groups, int *K, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *fdrmax, int *gapprox) {

/* Performs forward simulations for the gene classification problem using a Gamma/Gamma model */

/* Input arguments

   - B: number of forward simulations

   - Bsummary: number of simulations used to estimate the summary statistic. Less simulations are used if sufficient precision is achieved.

   - stepsum: summary statistic is computed up to 'stepsum' steps ahead in time

   - J: time horizon i.e. mpred*(J-Jini) obs per group are drawn from the predictive in each forward simulation

   - Jini: initial time

   - mpred: number of observations per groups to simulate from the predictive at each time point

   - util: util==1 means calling minfnrstfdr, util==2 calls maxwtpfp, util==3 calls maxamountde

   - cf: coefficients for the utility function. Ignored if util==1.

   - cf_sam: sampling cost. Ignored if util==1.

   - genelimit: only the genelimit genes with lowest v0 (prob of being equally expressed across all groups) will be used

   - v0thre: only genes with v0thre<v0 (prob of being equally expressed across all groups) will be used

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to (length ncol)

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - alpha0,nu: prior for lambda is Gamma(alpha0,nu)

   - balpha,nualpha: prior for alpha is Gamma(balpha,nualpha)

   - equalcv: set to 1 to indicate constant CV across groups, 0 otherwise

   - nclust: number of clusters in hyper-prior for (alpha,lambda)

   - cluslist: list of clusters to be taken into account when computing post prob. This option can be used to exclude clusters with very small probabilities

   - rho: mixing probabilities for clusters in hyper-prior for (alpha,lambda)

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - fdrmax: upper bound for E(FDR). Ignored if util!=1.

*/

/* Output arguments (each a vector of length B*(J-Jini) except for summary, a matrix (B*(J-Jini) rows, stepsum cols))

   - simid: simulation number

   - time: time

   - u: expected terminal utility (doesn't include sampling cost)

   - fdr: expected fdr

   - fnr: expected fnr

   - power: power as estimated by E(TP)/E(positives)

   - summary: summary statistic whose values will be used to define stopping regions. It's a matrix with the increment in utility up to 'stepsum' steps ahead in time.

*/



  int i,j,k, one=1, zero=0, ncolsumx, ncolxpred, uselpred, init0, *d, *dpred, *groupspred, deltat, nsel, *sel, *cluslist;

  double *sumx, *sumxpred, *sumxtot, *prodx, *prodxpred, *prodxtot, *nobsx, *nobsxpred, *nobsxtot, *xpred, *apred, *lpred, *m, *s, *v, *vpred, uobs, fdrobs, fnrobs, powobs, threshold, lhood, preceps=0.01;



  ncolxpred= (*mpred)*(*K);

  for (i=0, ncolsumx=0; i<(*npat); i++) { ncolsumx += ngrouppat[i]; }

  groupspred= ivector(0,ncolxpred);

  for (i=0; i<ncolxpred; i++) { groupspred[i]= floor(i/(*mpred)); }



  //Allocate memory

  cluslist= ivector(0,*nclust);

  sumx= dvector(0, (*nrow)*ncolsumx); sumxpred= dvector(0, (*nrow)*ncolsumx); sumxtot= dvector(0, (*nrow)*ncolsumx);

  if (*equalcv == 1) {

    prodx= dvector(0, (*nrow)); prodxpred= dvector(0, (*nrow)); prodxtot= dvector(0, (*nrow));

  } else {

    prodx= dvector(0, (*nrow)*ncolsumx); prodxpred= dvector(0, (*nrow)*ncolsumx); prodxtot= dvector(0, (*nrow)*ncolsumx);

  }

  nobsx= dvector(0, ncolsumx); nobsxpred= dvector(0, ncolsumx); nobsxtot= dvector(0, ncolsumx);

  xpred= dvector(0,(*nrow)*ncolxpred); dpred= ivector(0, *nrow); lpred= dvector(0, (*nrow)*(*K));

  if (*equalcv == 1) apred= dvector(0,*nrow); else apred= dvector(0,(*nrow)*(*K));

  v= dvector(0, (*nrow)*(*npat)); vpred= dvector(0, (*nrow)*(*npat));

  d= ivector(0, *nrow);

  m= dvector(0,*stepsum); s= dvector(0,*stepsum);



  //Compute sufficient statistics & Posterior probabilities of each expression pattern

  sel= ivector(0,*nrow);

  for (i=0; i<(*nrow); i++) { sel[i]= i; }

  compute_sumxC(sumx,prodx,nobsx,equalcv,nrow,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one);

  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;

  pp_ggC(v,&lhood,nrow,sel,ncol,x,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&zero,gapprox);



  //Select most DE genes

  sel_mostDEgenes(&nsel,sel,genelimit,v0thre,nrow,npat,v);



  //Loop over forward simulations

  for (i=0; i<(*B); i++) {

    uselpred= 0; init0= 1;

    copy_sumxC(sumxtot,prodxtot,nobsxtot,equalcv,&nsel,sel,&ncolsumx,sumx,prodx,nobsx); //copy contents of sumx, nobsx into sumxtot, nobsxtot



    //Loop over time

    for (j=0; j<((*J)-(*Jini)); j++) {

      simpred_ggC(xpred,dpred,apred,lpred,&uselpred,mpred,groups,K,&nsel,sel,nrow,ncol,x,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,v,npat,patterns,ngrouppat,sumx,prodx,nobsx,&one,gapprox);

      compute_sumxC(sumxpred,prodxpred,nobsxpred,equalcv,&nsel,sel,&ncolsumx,&ncolxpred,xpred,groupspred,K,npat,patterns,ngrouppat,&init0);

      compute_sumxC(sumxtot,prodxtot,nobsxtot,equalcv,&nsel,sel,&ncolsumx,&ncolxpred,xpred,groupspred,K,npat,patterns,ngrouppat,&zero);

      uselpred= 1; init0= 0;

      pp_ggC(vpred,&lhood,&nsel,sel,ncol,x,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox);

      //Optimal terminal decision and expected utility

      utgene_parC(&uobs, d, &fdrobs, &fnrobs, &powobs, &threshold, util, cf, &nsel, sel, vpred, npat, fdrmax);



      //Compute summary statistic

      if (j<(*J - *Jini -1)) {      //if we're not at the last time point compute summary stat

	for (k=0; k<(*stepsum); k++) {

	  deltat= (*mpred)*(k+1);

	  utgene_predC(m+k,s+k,&deltat,Bsummary,&preceps,util,cf,genelimit,v0thre,&nsel,sel,&one,nrow,ncol,x,groups,K,vpred,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,prob,npat,patterns,ngrouppat,fdrmax,sumxtot,prodxtot,nobsxtot,&one,gapprox);

          m[k]= m[k] - uobs - (k+1)*(*cf_sam)*((*util)!=1);   //for util!=1 substract sampling cost

	}

      } else {                      //if it's last time point don't compute summary (-9999 indicates missing)

	for (k=0; k<(*stepsum); k++) { m[k]= -9999; }

      }  



      //Observed utility if we were to stop now

      //uobsgeneC(&uobs, &fdrobs, &fnrobs, util, &nsel, sel, d, dpred, ltrue, cf);   



      //Save values

      simid[i*(*J - *Jini)+j]= i; time[i*(*J - *Jini)+j]= (*Jini)+j+1;

      u[i*(*J - *Jini)+j]= uobs;

      fdr[i*(*J - *Jini)+j]= fdrobs; fnr[i*(*J - *Jini)+j]= fnrobs; power[i*(*J - *Jini)+j]= powobs;

      for (k=0; k<(*stepsum); k++) { summary[(*stepsum)*(i*(*J - *Jini)+j)+k]= m[k]; }



    }  //End loop over time

  }  //End loop over forw simulations



  //Free memory

  free_ivector(groupspred,0,ncolxpred);

  free_ivector(cluslist,0,*nclust);

  free_dvector(sumx,0,(*nrow)*ncolsumx); free_dvector(sumxpred,0,(*nrow)*ncolsumx); free_dvector(sumxtot,0,(*nrow)*ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodx,0,(*nrow)); free_dvector(prodxpred,0,(*nrow)); free_dvector(prodxtot,0,(*nrow));

  } else {

    free_dvector(prodx,0,(*nrow)*ncolsumx); free_dvector(prodxpred,0,(*nrow)*ncolsumx); free_dvector(prodxtot,0,(*nrow)*ncolsumx);

  }

  free_dvector(nobsx,0,ncolsumx); free_dvector(nobsxpred,0,ncolsumx); free_dvector(nobsxtot,0,ncolsumx);

  free_dvector(xpred,0,(*nrow)*ncolxpred); free_ivector(dpred,0,*nrow); free_dvector(lpred,0,(*nrow)*(*K));

  if (*equalcv == 1) free_dvector(apred,0,*nrow); else free_dvector(apred,0,(*nrow)*(*K));



  free_dvector(v,0,(*nrow)*(*npat)); free_dvector(vpred,0,(*nrow)*(*npat));

  free_ivector(sel,0,*nrow);

  free_ivector(d,0, *nrow);

  free_dvector(m,0,*stepsum); free_dvector(s,0,*stepsum);



}





void forwsim_sampleC(int *simid, int *time, double *u, double *cc, double *summary, int *B, int *Bsummary, int *stepsum, int *J, int *Jini, int *mpred, int *genelimit, double *v0thre, int *nrow, int *ncol, double *x, int *groups, int *K, double *Kprob, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, int *gapprox) {

/* Performs forward simulations for the sample classification problem using a Gamma/Gamma model */

/* Input arguments

   - B: number of forward simulations

   - Bsummary: number of simulations used to estimate the summary statistic. Less simulations are used if sufficient precision is achieved.

   - stepsum: summary statistic is computed up to 'stepsum' steps ahead in time

   - J: time horizon i.e. mpred*(J-Jini) obs per group are drawn from the predictive in each forward simulation

   - Jini: initial time

   - mpred: number of observations per groups to simulate from the predictive at each time point

   - genelimit: only the genelimit genes with lowest v0 (prob of being equally expressed across all groups) will be used to classify samples

   - v0thre: only the genes with v0thre>=v0 (prob of being equally expressed across all groups) will be used to classify samples

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to (length ncol)

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - Kprob: vector with prior probabilities for each group

   - alpha0,nu: prior for lambda is Gamma(alpha0,nu)

   - balpha,nualpha: prior for alpha is Gamma(balpha,nualpha)

   - equalcv: set to 1 to indicate constant CV across groups, 0 otherwise

   - nclust: number of clusters in hyper-prior for (alpha,lambda)

   - cluslist: list of clusters to be taken into account when computing post prob. This option can be used to exclude clusters with very small probabilities

   - rho: mixing probabilities for clusters in hyper-prior for (alpha,lambda)

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used

*/

/* Output arguments

   - simid: simulation number (vector of length B*(J-Jini))

   - time: time (vector of length B*(J-Jini))

   - u: estimated overall prob. of correct classification (vector of length B*(J-Jini))

   - cc: matrix with estimated prob. of correct classification for each of the K groups (B*(J-Jini) rows, K cols)

   - summary: matrix with summary statistics whose values will be used to define stopping regions. It's a matrix with the increment in utility up to 'stepsum' steps ahead in time (B*(J-Jini) rows, stepsum cols)

*/



  int i,j,k, one=1, zero=0, ncolsumx, ncolxpred, uselpred, init0, *dpred, *groupspred, deltat, nsel=0, *sel, *ngroup, *nccpred, *cluslist;

  double *sumx, *sumxpred, *sumxtot, *prodx, *prodxpred, *prodxtot, *nobsx, *nobsxpred, *nobsxtot, *xpred, *apred, *lpred, *m, *s, *v, *vpred, uobs, *ccobs, *ccpred, lhood, preceps=0.001, seobs, sepred;



  ncolxpred= (*mpred)*(*K);

  for (i=0, ncolsumx=0; i<(*npat); i++) { ncolsumx += ngrouppat[i]; }

  groupspred= ivector(0,ncolxpred);

  for (i=0; i<ncolxpred; i++) { groupspred[i]= floor(i/(*mpred)); }



  //Allocate memory

  cluslist= ivector(0,*nclust);

  sumx= dvector(0, (*nrow)*ncolsumx); sumxpred= dvector(0, (*nrow)*ncolsumx); sumxtot= dvector(0, (*nrow)*ncolsumx);

  if (*equalcv == 1) {

    prodx= dvector(0, (*nrow)); prodxpred= dvector(0, (*nrow)); prodxtot= dvector(0, (*nrow));

  } else {

    prodx= dvector(0, (*nrow)*ncolsumx); prodxpred= dvector(0, (*nrow)*ncolsumx); prodxtot= dvector(0, (*nrow)*ncolsumx);

  }

  nobsx= dvector(0, ncolsumx); nobsxpred= dvector(0, ncolsumx); nobsxtot= dvector(0, ncolsumx);

  xpred= dvector(0,(*nrow)*ncolxpred); dpred= ivector(0, *nrow); lpred= dvector(0, (*nrow)*(*K));

  if (*equalcv == 1) apred= dvector(0,*nrow); else apred= dvector(0,(*nrow)*(*K));

  v= dvector(0, (*nrow)*(*npat)); vpred= dvector(0, (*nrow)*(*npat));

  sel= ivector(0,*nrow);

  ngroup= ivector(0,*K);



  ccobs= dvector(0,*K); ccpred= dvector(0,*K); nccpred= ivector(0,*K);

  m= dvector(0,*stepsum); s= dvector(0,*stepsum);



  //Compute sufficient statistics & Posterior probabilities of each expression pattern

  for (i=0; i<(*nrow); i++) { sel[i]= i; }

  compute_sumxC(sumx,prodx,nobsx,equalcv,nrow,sel,&ncolsumx, ncol, x, groups, K, npat, patterns, ngrouppat, &one);

  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;

  pp_ggC(v,&lhood,nrow,sel,ncol,x,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&zero,gapprox);



  //Select most DE genes

  sel_mostDEgenes(&nsel,sel,genelimit,v0thre,nrow,npat,v);



  //Loop over forward simulations

  for (i=0; i<(*B); i++) {

    uselpred= 0; init0= 1;

    copy_sumxC(sumxtot,prodxtot,nobsxtot,equalcv,&nsel,sel,&ncolsumx,sumx,prodx,nobsx); //copy contents of sumx, nobsx into sumxtot, nobsxtot



    //Loop over time

    for (j=0; j<((*J)-(*Jini)); j++) {

      simpred_ggC(xpred,dpred,apred,lpred,&uselpred,mpred,groups,K,&nsel,sel,nrow,ncol,x,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,v,npat,patterns,ngrouppat,sumx,prodx,nobsx,&one,gapprox);

      compute_sumxC(sumxpred,prodxpred,nobsxpred,equalcv,&nsel,sel,&ncolsumx,&ncolxpred,xpred,groupspred,K,npat,patterns,ngrouppat,&init0);

      compute_sumxC(sumxtot,prodxtot,nobsxtot,equalcv,&nsel,sel,&ncolsumx,&ncolxpred,xpred,groupspred,K,npat,patterns,ngrouppat,&zero);

      uselpred= 1; init0= 0;

      pp_ggC(vpred,&lhood,&nsel,sel,ncol,x,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox);



      //Prob. of correct classification if we were to stop collecting data at this time

      utsample_ggC(&uobs,&seobs,ccobs,ngroup,Bsummary,&preceps,genelimit,v0thre,&nsel,sel,&one,nrow,&ncolxpred,xpred,groupspred,vpred,K,Kprob,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,prob,npat,patterns,ngrouppat,&ncolsumx,sumxtot,prodxtot,nobsxtot,&one,gapprox);



      //Compute summary statistic

      if (j<(*J - *Jini -1)) {      //if we're not at the last time point compute summary stat

	for (k=0; k<(*stepsum); k++) {

	  deltat= (*mpred)*(k+1);

	  utsample_predC(m+k,&sepred,ccpred,nccpred,&deltat,Bsummary,&preceps,genelimit,v0thre,&nsel,sel,&one,nrow,ncol,x,groups,K,Kprob,vpred,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,prob,npat,patterns,ngrouppat,sumxtot,prodxtot,nobsxtot,&one,gapprox);

          m[k]= m[k] - uobs;

	}

      } else {                      //if it's last time point don't compute summary (-9999 indicates missing)

	for (k=0; k<(*stepsum); k++) { m[k]= -9999; }

      }  



      //Save values

      simid[i*(*J - *Jini)+j]= i; time[i*(*J - *Jini)+j]= (*Jini)+j+1;

      u[i*(*J - *Jini)+j]= uobs;

      for (k=0; k<(*K); k++) { cc[(*K)*(i*(*J - *Jini)+j)+k]= ccobs[k]; }

      for (k=0; k<(*stepsum); k++) { summary[(*stepsum)*(i*(*J - *Jini)+j)+k]= m[k]; }



    }  //End loop over time

  }  //End loop over forw simulations



  //Free memory

  free_ivector(groupspred,0,ncolxpred); free_ivector(cluslist,0,*nclust);

  free_dvector(sumx,0,(*nrow)*ncolsumx); free_dvector(sumxpred,0,(*nrow)*ncolsumx); free_dvector(sumxtot,0,(*nrow)*ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodx,0,(*nrow)); free_dvector(prodxpred,0,(*nrow)); free_dvector(prodxtot,0,(*nrow));

  } else {

    free_dvector(prodx,0,(*nrow)*ncolsumx); free_dvector(prodxpred,0,(*nrow)*ncolsumx); free_dvector(prodxtot,0,(*nrow)*ncolsumx);

  }

  free_dvector(nobsx,0,ncolsumx); free_dvector(nobsxpred,0,ncolsumx); free_dvector(nobsxtot,0,ncolsumx);

  free_dvector(xpred,0,(*nrow)*ncolxpred); free_ivector(dpred,0,*nrow); free_dvector(lpred,0,(*nrow)*(*K));

  if (*equalcv == 1) free_dvector(apred,0,*nrow); else free_dvector(apred,0,(*nrow)*(*K));



  free_dvector(v,0,(*nrow)*(*npat)); free_dvector(vpred,0,(*nrow)*(*npat));

  free_ivector(sel,0,*nrow);

  free_dvector(ccobs,0,*K); free_dvector(ccpred,0,*K); free_ivector(nccpred,0,*K);

  free_dvector(m,0,*stepsum); free_dvector(s,0,*stepsum);

  free_ivector(ngroup,0,*K);



}





/*********************************************************************************************

                                GENERALIZED GAMMA/GAMMA MODEL FIT

*********************************************************************************************/



void initpar_ggC(double *alpha0, double *nu, double *balpha, double *nualpha,  double *rho, double *prob, double *lhood, int *B, double *a_alpha0, double *b_alpha0, double *a_nu, double *b_nu, double *a_balpha, double *b_balpha, double *a_nualpha, double *b_nualpha, double *p_rho, double *p_prob, double *alpha0ini, double *nuini, double *balphaini, double *nualphaini, double *rhoini, double *probini, int *nrow, int *ncol, double *x, int *groups, int *K, int *equalcv, int *nclust, int *npat, int *patterns, int *ngrouppat, int *gapprox, int *trace) {

/* Refines initial parameter estimates for GaGa model by iteratively computing:

   (i) Posterior expectation of (balpha,nualpha,alpha0,nu) given probpat, probclus and the data

   (ii) Posterior expectation of probpat given (balpha,nualpha,alpha0,nu), probclus and the data

   (iii) Posterior expectatin of probclus given (balpha,nualpha,alpha0,nu), probpat and the data

   INPUT: as in fit_ggC. *B indicates number of MCMC iterations

*/

  int i, j, ncolsumx, one=1, *sel, *cluslist, settomean=1;

  double *sumx, *prodx, *nobsx, *sumxpred, *prodxpred, *nobsxpred, *sumd, *sumci, *sumalpha, *sumlogalpha, *suminvlambda, *sumlambda, *sumloglambda, *v, *ngroupstot;



  cluslist= ivector(0,*nclust);

  sumalpha= dvector(0,*nclust); sumlogalpha= dvector(0,*nclust);

  suminvlambda= dvector(0,*nclust); sumlambda= dvector(0,*nclust); sumloglambda= dvector(0,*nclust);

  ngroupstot= dvector(0,*nclust);



  //Pre-compute sufficient statistics

  sel= ivector(0,*nrow);

  for (i=0; i<(*nrow); i++) { sel[i]= i; }

  for (i=0, ncolsumx=0; i<*npat; i++) { ncolsumx += ngrouppat[i]; }

  sumx= dvector(0,(*nrow)*ncolsumx); nobsx= dvector(0,ncolsumx);

  sumxpred= dvector(0,(*nrow)*ncolsumx); nobsxpred= dvector(0,ncolsumx);

  if (*equalcv == 1) {

    prodx= dvector(0,(*nrow));

    prodxpred= dvector(0,(*nrow));

    for (i=0; i<(*nrow)*ncolsumx; i++) { sumxpred[i]= 0; }

    for (i=0; i<(*nrow); i++) { prodxpred[i]= 0; }

  } else {

    prodx= dvector(0,(*nrow)*ncolsumx);

    prodxpred= dvector(0,(*nrow)*ncolsumx);

    for (i=0; i<(*nrow)*ncolsumx; i++) { sumxpred[i]= prodxpred[i]= 0; }

  }

  for (i=0; i<ncolsumx; i++) { nobsxpred[i]= 0; }

  compute_sumxC(sumx,prodx,nobsx,equalcv,nrow,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one);



  //MCMC

  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;

  v= dvector(0,(*nrow)*(*npat));

  sumd= dvector(0,*npat); sumci= dvector(0,*nclust);

  balpha[0]= balphaini[0]; nualpha[0]= nualphaini[0];

  for (i=0; i<*nclust; i++) { alpha0[i]= alpha0ini[i]; nu[i]= nuini[i]; }

  for (i=0; i<*npat; i++) { prob[i]= probini[i]; }

  for (i=0; i<*nclust; i++) { rho[i]= rhoini[i]; }

  for (i=1; (i<*B); i++) {

    pp_ggC(v,lhood+i-1,nrow,sel,ncol,x,groups,K,alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,equalcv,nclust,cluslist,rho+(i-1)*(*nclust),prob+(i-1)*(*npat),npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox);

    simpar_ggC(ngroupstot,sumd,sumci,sumalpha,sumlogalpha,suminvlambda,sumlambda,sumloglambda,groups,K,nrow,alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,equalcv,nclust,rho+(i-1)*(*nclust),v,npat,patterns,ngrouppat,sumx,prodx,nobsx,gapprox);

    simhyperpar_ggC(alpha0+i*(*nclust),nu+i*(*nclust),balpha+i,nualpha+i,nclust,rho+i*(*nclust),prob+i*(*npat),cluslist,equalcv,a_alpha0,b_alpha0,a_nu,b_nu,a_balpha,b_balpha,a_nualpha,b_nualpha,p_rho,p_prob,nrow,sumd,sumci,ngroupstot,sumalpha,sumlogalpha,suminvlambda,sumlambda,sumloglambda,npat,ngrouppat,gapprox,&settomean);



    if ((*trace==1) && (i>1)) {

      for (j=0; j<*nclust; j++) printf("%.6f ",alpha0[(i-1)*(*nclust)+j]);

      for (j=0; j<*nclust; j++) printf("%.6f ",nu[(i-1)*(*nclust)+j]);

      printf("%.6f %.6f ",balpha[0],nualpha[0]);

      for (j=0; j<*npat; j++) printf("%.6f ",prob[(i-1)*(*npat)+j]);

      for (j=0; j<*nclust; j++) printf("%.6f ",rho[(i-1)*(*nclust)+j]);

      printf("%.6f",lhood[i-1]);

      printf("\n");

    }

  }

  pp_ggC(v,lhood+i-1,nrow,sel,ncol,x,groups,K,alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,equalcv,nclust,cluslist,rho+(i-1)*(*nclust),prob+(i-1)*(*npat),npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox); //compute log-likelihood at last iteration



  for (j=0; j<*nclust; j++) printf("%.6f ",alpha0[(i-1)*(*nclust)+j]);

  for (j=0; j<*nclust; j++) printf("%.6f ",nu[(i-1)*(*nclust)+j]);

  printf("%.6f %.6f ",balpha[0],nualpha[0]);

  for (j=0; j<*npat; j++) printf("%.6f ",prob[(i-1)*(*npat)+j]);

  for (j=0; j<*nclust; j++) printf("%.6f ",rho[(i-1)*(*nclust)+j]);

  printf("%.6f \n",lhood[i-1]);



  free_ivector(sel,0,*nrow);

  free_dvector(sumd,0,*npat); free_dvector(sumci,0,*nclust); free_dvector(v,0,(*nrow)*(*npat));

  free_dvector(sumx,0,(*nrow)*ncolsumx); free_dvector(nobsx,0,ncolsumx);

  free_dvector(sumxpred,0,(*nrow)*ncolsumx);  free_dvector(nobsxpred,0,ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodx,0,(*nrow));

    free_dvector(prodxpred,0,(*nrow));

  } else {

    free_dvector(prodx,0,(*nrow)*ncolsumx);

    free_dvector(prodxpred,0,(*nrow)*ncolsumx);

  }



  free_ivector(cluslist,0,*nclust);

  free_dvector(sumalpha,0,*nclust); free_dvector(sumlogalpha,0,*nclust);

  free_dvector(suminvlambda,0,*nclust); free_dvector(sumlambda,0,*nclust); 

  free_dvector(sumloglambda,0,*nclust); free_dvector(ngroupstot,0,*nclust);

}



void fit_ggC(double *alpha0, double *nu, double *balpha, double *nualpha,  double *rho, double *prob, double *lhood, int *B, double *a_alpha0, double *b_alpha0, double *a_nu, double *b_nu, double *a_balpha, double *b_balpha, double *a_nualpha, double *b_nualpha, double *p_rho, double *p_prob, double *alpha0ini, double *nuini, double *balphaini, double *nualphaini, double *rhoini, double *probini, int *nrow, int *ncol, double *x, int *groups, int *K, int *equalcv, int *nclust, int *npat, int *patterns, int *ngrouppat, int *gapprox, int *trace) {

/* Fits GaGa model via Gibbs' sampling */

/* Input

   - B: number of MCMC samples to draw

   - b_alpha0: prior for alpha0 is Gamma(a_alpha0,b_alpha0)

   - a_nu, b_nu: prior for nu is Gamma(a_nu,b_nu)

   - a_balpha, b_balpha: prior for balpha is Gamma(a_balpha,b_balpha)

   - a_nualpha, b_nualpha: prior for nualpha is Gamma(a_nualpha,b_nualpha)

   - p_rho: prior for rho is Dirichlet(p_rho). (vector of length nrho)

   - p_prob: prior for prob is Dirichlet(p_prob). (vector of length npat)

   - alpha0ini: vector of length nclust with initial values for hyper-parameter alpha0

   - nuini: vector of length nclust with initial values for hyper-parameter nu

   - balphaini: vector of length nclust with initial values for hyper-parameter balpha

   - nualphaini: vector of length nclust with initial values for hyper-parameter nualpha

   - rhoini: initial value for hyper-parameter rho

   - probini: initial value for hyper-parameter prob

   - nrow: number of rows of x (genes). 

   - ncol: number of columns of x (samples).

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - equalcv: set to 1 to indicate constant CV across groups, 0 otherwise

   - nclust:

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used

   - trace: for trace==1 iteration progress is printed

   Output

   - alpha0: matrix (*B rows, *nclust cols) with posterior draws for hyper-parameter alpha0

   - nu: matrix (*B rows, *nclust cols) with posterior draws for hyper-parameter nu

   - balpha: vector (length *B) with posterior draws for hyper-parameter balpha

   - nualpha: vector (length *B) with posterior draws for hyper-parameter nualpha

   - rho: matrix (*B rows, *nclust cols) with posterior draws for hyper-parameter rho

   - prob: matrix with *B rows and *npat columns with posterior draws for hyper-parameter prob

   - lhood: vector (length *B) with log-likelihood evaluated at posterior draws

*/



  int i, ncolsumx, one=1, *sel, B10, *cluslist, settomean=0;

  double *sumx, *prodx, *nobsx, *sumxpred, *prodxpred, *nobsxpred, *sumd, *sumci, *sumalpha, *sumlogalpha, *suminvlambda, *sumlambda, *sumloglambda, *v, *ngroupstot;



  if (*B < 10) (*B)= 10;  //ensure a minimum amount of iterations

  cluslist= ivector(0,*nclust);

  sumalpha= dvector(0,*nclust); sumlogalpha= dvector(0,*nclust);

  suminvlambda= dvector(0,*nclust); sumlambda= dvector(0,*nclust); sumloglambda= dvector(0,*nclust);

  ngroupstot= dvector(0,*nclust);



  //Pre-compute sufficient statistics

  sel= ivector(0,*nrow);

  for (i=0; i<(*nrow); i++) { sel[i]= i; }

  for (i=0, ncolsumx=0; i<*npat; i++) { ncolsumx += ngrouppat[i]; }

  sumx= dvector(0,(*nrow)*ncolsumx); nobsx= dvector(0,ncolsumx);

  sumxpred= dvector(0,(*nrow)*ncolsumx); nobsxpred= dvector(0,ncolsumx);

  if (*equalcv == 1) {

    prodx= dvector(0,(*nrow));

    prodxpred= dvector(0,(*nrow));

    for (i=0; i<(*nrow)*ncolsumx; i++) { sumxpred[i]= 0; }

    for (i=0; i<(*nrow); i++) { prodxpred[i]= 0; }

  } else {

    prodx= dvector(0,(*nrow)*ncolsumx);

    prodxpred= dvector(0,(*nrow)*ncolsumx);

    for (i=0; i<(*nrow)*ncolsumx; i++) { sumxpred[i]= prodxpred[i]= 0; }

  }

  for (i=0; i<ncolsumx; i++) { nobsxpred[i]= 0; }

  compute_sumxC(sumx,prodx,nobsx,equalcv,nrow,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one);



  //MCMC

  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;

  v= dvector(0,(*nrow)*(*npat));

  sumd= dvector(0,*npat); sumci= dvector(0,*nclust);

  balpha[0]= balphaini[0]; nualpha[0]= nualphaini[0];

  for (i=0; i<*nclust; i++) { alpha0[i]= alpha0ini[i]; nu[i]= nuini[i]; }

  for (i=0; i<*npat; i++) { prob[i]= probini[i]; }

  for (i=0; i<*nclust; i++) { rho[i]= rhoini[i]; }

  if ((*trace) ==1) printf("Starting Gibbs sampler... \n");

  for (i=1, B10=(*B)/10; i<*B; i++) {

    pp_ggC(v,lhood+i-1,nrow,sel,ncol,x,groups,K,alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,equalcv,nclust,cluslist,rho+(i-1)*(*nclust),prob+(i-1)*(*npat),npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox);

    simpar_ggC(ngroupstot,sumd,sumci,sumalpha,sumlogalpha,suminvlambda,sumlambda,sumloglambda,groups,K,nrow,alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,equalcv,nclust,rho+(i-1)*(*nclust),v,npat,patterns,ngrouppat,sumx,prodx,nobsx,gapprox);

    simhyperpar_ggC(alpha0+i*(*nclust),nu+i*(*nclust),balpha+i,nualpha+i,nclust,rho+i*(*nclust),prob+i*(*npat),cluslist,equalcv,a_alpha0,b_alpha0,a_nu,b_nu,a_balpha,b_balpha,a_nualpha,b_nualpha,p_rho,p_prob,nrow,sumd,sumci,ngroupstot,sumalpha,sumlogalpha,suminvlambda,sumlambda,sumloglambda,npat,ngrouppat,gapprox,&settomean);

    if ((*trace==1) && ((i % B10)==0)) printf("  %d iterations \n",i);

  }

  pp_ggC(v,lhood+i-1,nrow,sel,ncol,x,groups,K,alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,equalcv,nclust,cluslist,rho+(i-1)*(*nclust),prob+(i-1)*(*npat),npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox); //compute log-likelihood at last iteration



  free_ivector(sel,0,*nrow);

  free_dvector(sumd,0,*npat); free_dvector(sumci,0,*nclust); free_dvector(v,0,(*nrow)*(*npat));

  free_dvector(sumx,0,(*nrow)*ncolsumx); free_dvector(nobsx,0,ncolsumx);

  free_dvector(sumxpred,0,(*nrow)*ncolsumx);  free_dvector(nobsxpred,0,ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodx,0,(*nrow));

    free_dvector(prodxpred,0,(*nrow));

  } else {

    free_dvector(prodx,0,(*nrow)*ncolsumx);

    free_dvector(prodxpred,0,(*nrow)*ncolsumx);

  }



  free_ivector(cluslist,0,*nclust);

  free_dvector(sumalpha,0,*nclust); free_dvector(sumlogalpha,0,*nclust);

  free_dvector(suminvlambda,0,*nclust); free_dvector(sumlambda,0,*nclust); 

  free_dvector(sumloglambda,0,*nclust); free_dvector(ngroupstot,0,*nclust);

}



double fnopat(double *th) {

  double lhood;

  lhoodnopat(esteppars.probest,esteppars.rhoest,&lhood,esteppars.justlhood,esteppars.nsel,esteppars.sel,esteppars.K,th+1,th+ *esteppars.nclust+1,th+2*(*esteppars.nclust)+1,th+2*(*esteppars.nclust)+2,esteppars.equalcv,esteppars.nclust,esteppars.rho,esteppars.prob,esteppars.npat,esteppars.patterns,esteppars.ngrouppat,esteppars.colini,esteppars.ncolsumx,esteppars.sumx,esteppars.prodx,esteppars.nobsx,esteppars.gapprox);

    return(-lhood);

}

double f(double *th) { 

  int i;

  double lhood;

  //  __iglobal++;  //debug

  //printf("%d \n",__iglobal); //debug



  //avoid numerical overflow by truncating pars at exp(700)  

  for(i=1; i<= (*esteppars.nclust); i++) { 

    if (th[i]>700) th[i]= 700.0;                                           //alpha0

    if (th[i + *esteppars.nclust]>700) th[i + *esteppars.nclust]= 700.0;   //nu

  }

  if (th[2*(*esteppars.nclust)+1]>700) th[2*(*esteppars.nclust)+1]= 700.0; //balpha

  if (th[2*(*esteppars.nclust)+2]>700) th[2*(*esteppars.nclust)+2]= 700.0; //nualpha



  estep(esteppars.probest,esteppars.rhoest,&lhood,esteppars.justlhood,esteppars.nsel,esteppars.sel,esteppars.K,th+1,th+ *esteppars.nclust+1,th+2*(*esteppars.nclust)+1,th+2*(*esteppars.nclust)+2,esteppars.equalcv,esteppars.nclust,esteppars.rho,esteppars.prob,esteppars.npat,esteppars.patterns,esteppars.ngrouppat,esteppars.colini,esteppars.ncolsumx,esteppars.sumx,esteppars.prodx,esteppars.nobsx,esteppars.gapprox);

  return(-lhood);

}



void fitEM_ggC(double *alpha0, double *nu, double *balpha, double *nualpha,  double *rho, double *prob, double *lhood, int *B, double *alpha0ini, double *nuini, double *balphaini, double *nualphaini, double *rhoini, double *probini, int *nrow, int *ncol, double *x, int *groups, int *K, int *equalcv, int *nclust, int *npat, int *patterns, int *ngrouppat, int *gapprox, int *trace)

 {

//Fit GaGa or MiGaGa model via EM algorithm

/* Input: as in fit_ggC, removing prior parameters a_alpha0,...,p_prob.

          B denotes the maximum number of iterations

   Output:

   - alpha0, nu, balpha, nualpha, rho, prob: parameter estimates

   - lhood: expected log-likelihood at the parameter estimate

*/



  int i, j, npar, iter, maxit, ncolsumx, *sel, one=1, justlhood, *colini;

  double *sumx, *prodx, *nobsx, *th, *thold, **dirth, ftol=1.0e-2, fret, fold, eps, epsth;





  //Pre-compute sufficient statistics

  colini= ivector(0,*npat);

  colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

  ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];



  sel= ivector(0,*nrow);

  for (i=0; i<(*nrow); i++) { sel[i]= i; }

  sumx= dvector(0,(*nrow)*ncolsumx); nobsx= dvector(0,ncolsumx);

  if (*equalcv == 1) {

    prodx= dvector(0,(*nrow));

  } else {

    prodx= dvector(0,(*nrow)*ncolsumx);

  }

  compute_sumxC(sumx,prodx,nobsx,equalcv,nrow,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one);



  //Initialize global variable esteppars with params that must be passed to estep

  esteppars.probest= dvector(0,*npat-1);

  esteppars.rhoest= dvector(0,*nclust-1);

  esteppars.lhood= dvector(0,1);

  esteppars.justlhood= &justlhood; 

  esteppars.nsel= nrow;

  esteppars.sel= sel;

  esteppars.K= K;

  esteppars.equalcv= equalcv;

  esteppars.nclust= nclust;

  esteppars.rho= rho;

  esteppars.prob= prob;

  esteppars.npat= npat;

  esteppars.patterns= patterns;

  esteppars.ngrouppat= ngrouppat;

  esteppars.colini= colini;

  esteppars.ncolsumx= &ncolsumx;

  esteppars.sumx= sumx;

  esteppars.prodx= prodx;

  esteppars.nobsx= nobsx;

  esteppars.gapprox= gapprox;



  //initial parameter values and initial directions for miminization algorithm

  npar= 2*(*nclust)+2;

  th= dvector(1,npar);

  thold= dvector(1,npar);

  for (i=1; i<= *nclust; i++) {

    th[i]= log(alpha0ini[i-1]);

    th[*nclust +i]= log(nuini[i-1]);

  }

  th[2*(*nclust)+1]= log(*balphaini);

  th[2*(*nclust)+2]= log(*nualphaini);



  for (i=0; i< *npat; i++) { esteppars.probest[i]= esteppars.prob[i]= probini[i]; }

  for (i=0; i< *nclust; i++) { esteppars.rhoest[i]= esteppars.rho[i]= rhoini[i]; }



  dirth= dmatrix(1,npar,1,npar);

  for (i=1; i<= npar; i++) {

    for (j=1; j<= npar; j++) {

      if (i==j) dirth[i][j]= 1; else dirth[i][j]= 0;

    }

  }



  //Refine initial estimates

  if (*trace ==1) printf("Refining initial estimates...\n");

  maxit= 5; 

  for (i=0, eps=1, justlhood=0; eps>.01; i++) {

    estep(prob,rho,lhood,&justlhood,nrow,sel,K,th+1,th+ *nclust +1,th+2*(*nclust)+1,th+2*(*nclust)+2,equalcv,nclust,esteppars.rho,esteppars.prob,npat,patterns,ngrouppat,colini,&ncolsumx,sumx,prodx,nobsx,gapprox);

    eps= 0;

    for (j=0; j< *npat; j++) { eps+= fabs(esteppars.probest[j]-prob[j]); esteppars.probest[j]= prob[j]; }

    for (j=0; j< *nclust; j++) { eps+= fabs(esteppars.rhoest[j]-rho[j]); esteppars.rhoest[j]= rho[j]; }

  }



  //EM algorithm

  if (*trace ==1) {

    printf("Starting EM algorithm...\n"); 

    for (j=1; j<= *nclust; j++) printf("alpha0[%d] ",j);

    for (j=1; j<= *nclust; j++) printf("    nu[%d] ",j);

    printf("   balpha    nualpha ");

    for (j=0; j< *npat; j++) printf("probpat[%d] ",j+1);

    if (*nclust>1) for (j=0; j< *nclust; j++) printf("probclus[%d] ",j+1);

    printf("logl\n");

    for (j=1; j<= *nclust; j++) printf("%9f ",exp(th[j]));

    for (j=1; j<= *nclust; j++) printf("%9f ",exp(th[*nclust +j]));

    printf("%9f %9f ",exp(th[2*(*nclust)+1]),exp(th[2*(*nclust)+2]));

    for (j=0; j< *npat; j++) printf("%10f ",prob[j]);

    if (*nclust>1) for (j=0; j< *nclust; j++) printf("%11f ",rho[j]);

    printf("%9f\n",*lhood);

  }



  fold= *lhood;

  for (i=0,eps=1,epsth=1; (i< *B) && (eps>.01) && (epsth>.01); i++) {

    justlhood= 1;

    for (j=1; j<(2*(*nclust)+2); j++) { thold[j]= th[j]; }

    minimize(th,dirth,npar,ftol,&iter,&fret,f,maxit);  //M-step



    justlhood= 0;

    estep(prob,rho,lhood,&justlhood,nrow,sel,K,th+1,th+ *nclust +1,th+2*(*nclust)+1,th+2*(*nclust)+2,equalcv,nclust,esteppars.rho,esteppars.prob,npat,patterns,ngrouppat,colini,&ncolsumx,sumx,prodx,nobsx,gapprox); //E-step

    eps= fabs(*lhood - fold);

    fold= *lhood;

    for (j=1, epsth=0; j<(2*(*nclust)+2); j++) { epsth+= (th[j]-thold[j])*(th[j]-thold[j]); }



    for (j=0; j< *npat; j++) { esteppars.probest[j]= prob[j]; }

    for (j=0; j< *nclust; j++) { esteppars.rhoest[j]= rho[j]; }



    if (*trace ==1) { 

      for (j=1; j<= *nclust; j++) printf("%9f ",exp(th[j]));

      for (j=1; j<= *nclust; j++) printf("%9f ",exp(th[*nclust +j]));

      printf("%9f %9f ",exp(th[2*(*nclust)+1]),exp(th[2*(*nclust)+2]));

      for (j=0; j< *npat; j++) printf("%10f ",prob[j]);

      if (*nclust>1) for (j=0; j< *nclust; j++) printf("%11f ",rho[j]);

      printf("%9f\n",*lhood);

    }

  }



  //Refine mixing probabilities after other hyper-pars have converged

  for (i=0, eps=1, justlhood=0; eps>.01; i++) {

    estep(prob,rho,lhood,&justlhood,nrow,sel,K,th+1,th+ *nclust +1,th+2*(*nclust)+1,th+2*(*nclust)+2,equalcv,nclust,esteppars.rho,esteppars.prob,npat,patterns,ngrouppat,colini,&ncolsumx,sumx,prodx,nobsx,gapprox);

    eps= 0;

    for (j=0; j< *npat; j++) { eps+= fabs(esteppars.probest[j]-prob[j]); esteppars.probest[j]= prob[j]; }

    for (j=0; j< *nclust; j++) { eps+= fabs(esteppars.rhoest[j]-rho[j]); esteppars.rhoest[j]= rho[j]; }

    if (*trace ==1) { 

      for (j=1; j<= *nclust; j++) printf("%9f ",exp(th[j]));

      for (j=1; j<= *nclust; j++) printf("%9f ",exp(th[*nclust +j]));

      printf("%9f %9f ",exp(th[2*(*nclust)+1]),exp(th[2*(*nclust)+2]));

      for (j=0; j< *npat; j++) printf("%10f ",prob[j]);

      if (*nclust>1) for (j=0; j< *nclust; j++) printf("%11f ",rho[j]);

      printf("%9f\n",*lhood);

    }

  }



  //save estimates in output arguments

  for (i=1; i<= *nclust; i++) { alpha0[i-1]= exp(th[i]); nu[i-1]= exp(th[*nclust +i]); }

  (*balpha)= exp(th[2*(*nclust)+1]); (*nualpha)= exp(th[2*(*nclust)+2]);



  free_dvector(esteppars.probest,0,*npat-1);

  free_dvector(esteppars.rhoest,0,*nclust-1);

  free_dvector(esteppars.lhood,0,1);

  free_dvector(th,1,npar); free_dvector(thold,1,npar);

  free_dmatrix(dirth,1,npar,1,npar);



  free_ivector(sel,0,*nrow);

  free_ivector(colini,0,*npat);

  free_dvector(sumx,0,(*nrow)*ncolsumx); free_dvector(nobsx,0,ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodx,0,(*nrow));

  } else {

    free_dvector(prodx,0,(*nrow)*ncolsumx);

  }



}



void fitMH_ggC(double *acprop, double *alpha0, double *nu, double *balpha, double *nualpha,  double *rho, double *prob, double *lhood, int *B, double *a_alpha0, double *b_alpha0, double *a_nu, double *b_nu, double *a_balpha, double *b_balpha, double *a_nualpha, double *b_nualpha, double *p_rho, double *p_prob, double *alpha0ini, double *nuini, double *balphaini, double *nualphaini, double *rhoini, double *probini, int *nrow, int *ncol, double *x, int *groups, int *K, int *equalcv, int *nclust, int *npat, int *patterns, int *ngrouppat, int *gapprox, int *trace, int *method, int *Bgibbs, double *h_alpha0, double *h_nu, double *h_balpha, double *h_nualpha, double *h_rho, double *h_prob) {

/* Fits GaGa model via Metropolis-Hastings or Simulated Annealing */

/* Input: as in fit_ggC, adding

   - method: 1 for MH, 0 for Simulated Annealing

   - Bgibbs: number of Gibb's sampling iterations to train the MH proposal parameters.

   - h_alpha0, h_nu, h_balpha, h_nualpha, h_rho, h_prob: MH proposal parameters. Ignored unless Bgibbs==0. See dproposal for the parameter meaning

   Output: as in fit_ggC, plus acprop: proportion of accepted proposals. Note: here lhood returns log(likelihood)+log(prior), in fit_ggC it just returns log(likelihood)

*/



  int i, j, ncolsumx, one=1, *sel, B10, *cluslist, logscale=1, Bburnin, gibbstrace=0;

  double *sumx, *prodx, *nobsx, *sumxpred, *prodxpred, *nobsxpred, *sumd, *sumci, *v, *ngroupstot, *a0gibbs, *nugibbs, *balphagibbs, *nualphagibbs, *rhogibbs, *probgibbs, *lgibbs, *vectemp, *vectemp2, propdens0, propdens1, lambda, u, temperature, maxlhood;



  if (*Bgibbs > 0) {      //Train MH proposal parameters

    if ((*trace)==1) printf("Training proposal parameters...\n");

    //Do a few Gibbs' iterations 

    a0gibbs= dvector(0,(*Bgibbs)*(*nclust)); nugibbs= dvector(0,(*Bgibbs)*(*nclust));

    balphagibbs= dvector(0,*Bgibbs); nualphagibbs= dvector(0,*Bgibbs);

    rhogibbs= dvector(0,(*Bgibbs)*(*nclust)); probgibbs= dvector(0,(*Bgibbs)*(*npat));

    lgibbs= dvector(0,*Bgibbs);

    fit_ggC(a0gibbs,nugibbs,balphagibbs,nualphagibbs,rhogibbs,probgibbs,lgibbs,Bgibbs,a_alpha0,b_alpha0,a_nu,b_nu,a_balpha,b_balpha,a_nualpha,b_nualpha,p_rho,p_prob,alpha0ini,nuini,balphaini,nualphaini,rhoini,probini,nrow,ncol,x,groups,K,equalcv,nclust,npat,patterns,ngrouppat,gapprox,&gibbstrace);



    //Set proposal parameters

    Bburnin= 0.2*(*Bgibbs);   //Gibbs samples to exclude from training sample

    vectemp= dvector(0,*nclust); 

    vectemp2= dvector(0, *nclust);

    colCV(vectemp,a0gibbs+Bburnin*(*nclust),*Bgibbs-Bburnin,*nclust);

    (*h_alpha0)= 1.0/pow(meanx(vectemp,*nclust -1),2);

    colCVinv(vectemp,nugibbs+Bburnin*(*nclust),*Bgibbs-Bburnin,*nclust); 

    (*h_nu)= 1.0/pow(meanx(vectemp,*nclust -1),2);



    (*h_balpha)= 1.0/pow(cv(balphagibbs,Bburnin,*Bgibbs -1),2);

    (*h_nualpha)= 1.0/pow(cvinv(nualphagibbs,Bburnin,*Bgibbs -1),2);

    free_dvector(vectemp,0,*nclust); free_dvector(vectemp2,0,*nclust);



    if (*nclust >1) {

      vectemp= dvector(0,*nclust); vectemp2= dvector(0, *nclust);

      colMeans(vectemp,rhogibbs+Bburnin*(*nclust),*Bgibbs-Bburnin,*nclust);

      colVar(vectemp2,rhogibbs+Bburnin*(*nclust),*Bgibbs-Bburnin,*nclust);

      for (i=0; i<*nclust; i++) { vectemp[i]= vectemp[i]*(1.0-vectemp[i])/vectemp2[i] - 1.0; }

      (*h_rho)= meanx(vectemp,*nclust -1);

      free_dvector(vectemp,0,*nclust); free_dvector(vectemp2,0,*nclust);

    } else { (*h_rho)= 1; }



    vectemp= dvector(0,*npat); vectemp2= dvector(0, *npat);

    colMeans(vectemp,probgibbs+Bburnin*(*npat),*Bgibbs - Bburnin,*npat);

    colVar(vectemp2,probgibbs+Bburnin*(*npat),*Bgibbs - Bburnin,*npat);

    for (i=0; i<*npat; i++) { vectemp[i]= vectemp[i]*(1.0-vectemp[i])/vectemp2[i] - 1.0; }

    (*h_prob)= meanx(vectemp,*npat -1);

    free_dvector(vectemp,0,*npat); free_dvector(vectemp2,0,*npat);



    free_dvector(a0gibbs,0,(*Bgibbs)*(*nclust)); free_dvector(nugibbs,0,(*Bgibbs)*(*nclust));

    free_dvector(balphagibbs,0,*Bgibbs); free_dvector(nualphagibbs,0,*Bgibbs);

    free_dvector(rhogibbs,0,(*Bgibbs)*(*nclust)); free_dvector(probgibbs,0,(*Bgibbs)*(*npat));

    free_dvector(lgibbs,0,*Bgibbs);



    balpha[0]= balphagibbs[*Bgibbs-1]; nualpha[0]= nualphagibbs[*Bgibbs-1];

    for (i=0; i<*nclust; i++) { alpha0[i]= a0gibbs[(*Bgibbs-1)*(*nclust)+i]; nu[i]= nugibbs[(*Bgibbs-1)*(*nclust)+i]; }

    for (i=0; i<*npat; i++) { prob[i]= probgibbs[(*Bgibbs-1)*(*npat)+i]; }

    for (i=0; i<*nclust; i++) { rho[i]= rhogibbs[(*Bgibbs-1)*(*nclust)+i]; }



    //(*h_alpha0)= .5*(*h_alpha0); (*h_nu)= .5*(*h_nu); (*h_balpha)= .5*(*h_balpha); (*h_nualpha)= .5*(*h_nualpha);

    //(*h_prob)= 2.0*(*h_prob); (*h_rho)= 2.0*(*h_rho);

    if ((*trace)==1) printf("Done.\n");

  } else {

    balpha[0]= balphaini[0]; nualpha[0]= nualphaini[0];

    for (i=0; i<*nclust; i++) { alpha0[i]= alpha0ini[i]; nu[i]= nuini[i]; }

    for (i=0; i<*npat; i++) { prob[i]= probini[i]; }

    for (i=0; i<*nclust; i++) { rho[i]= rhoini[i]; }

  }



  if (*B < 10) (*B)= 10;  //ensure a minimum amount of iterations

  cluslist= ivector(0,*nclust);

  ngroupstot= dvector(0,*nclust);



  //Pre-compute sufficient statistics

  sel= ivector(0,*nrow);

  for (i=0; i<(*nrow); i++) { sel[i]= i; }

  for (i=0, ncolsumx=0; i<*npat; i++) { ncolsumx += ngrouppat[i]; }

  sumx= dvector(0,(*nrow)*ncolsumx); nobsx= dvector(0,ncolsumx);

  sumxpred= dvector(0,(*nrow)*ncolsumx); nobsxpred= dvector(0,ncolsumx);

  if (*equalcv == 1) {

    prodx= dvector(0,(*nrow));

    prodxpred= dvector(0,(*nrow));

    for (i=0; i<(*nrow)*ncolsumx; i++) { sumxpred[i]= 0; }

    for (i=0; i<(*nrow); i++) { prodxpred[i]= 0; }

  } else {

    prodx= dvector(0,(*nrow)*ncolsumx);

    prodxpred= dvector(0,(*nrow)*ncolsumx);

    for (i=0; i<(*nrow)*ncolsumx; i++) { sumxpred[i]= prodxpred[i]= 0; }

  }

  for (i=0; i<ncolsumx; i++) { nobsxpred[i]= 0; }

  compute_sumxC(sumx,prodx,nobsx,equalcv,nrow,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one);



  //MCMC

  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;

  v= dvector(0,(*nrow)*(*npat));

  sumd= dvector(0,*npat); sumci= dvector(0,*nclust);

  pp_ggC(v,lhood,nrow,sel,ncol,x,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox);

  lhood[0]+= dprior_hyperpar(alpha0,nu,balpha,nualpha,rho,prob,nclust,npat,a_alpha0,b_alpha0,a_nu,b_nu,a_balpha,b_balpha,a_nualpha,b_nualpha,p_rho,p_prob,&logscale);



  if ((*trace) ==1) {

    if ((*method)==1) printf("Starting Metropolis-Hastings sampler... \n"); else printf("Starting Simulated Annealing. Printing best log-posterior found so far...\n  Iteration log-posterior\n");

  }

  for (i=1, B10=(*B)/10, *acprop= 0, maxlhood= lhood[0]; i<*B; i++) {

    rproposal(alpha0+i*(*nclust),nu+i*(*nclust),balpha+i,nualpha+i,rho+i*(*nclust),prob+i*(*npat),alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,rho+(i-1)*(*nclust),prob+(i-1)*(*npat),nclust,npat,h_alpha0,h_nu,h_balpha,h_nualpha,h_rho,h_prob); //generate proposal



    pp_ggC(v,lhood+i,nrow,sel,ncol,x,groups,K,alpha0+i*(*nclust),nu+i*(*nclust),balpha+i,nualpha+i,equalcv,nclust,cluslist,rho+i*(*nclust),prob+i*(*npat),npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox);

    lhood[i]+= dprior_hyperpar(alpha0+i*(*nclust),nu+i*(*nclust),balpha+i,nualpha+i,rho+i*(*nclust),prob+i*(*npat),nclust,npat,a_alpha0,b_alpha0,a_nu,b_nu,a_balpha,b_balpha,a_nualpha,b_nualpha,p_rho,p_prob,&logscale);

    propdens0= dproposal(alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,rho+(i-1)*(*nclust),prob+(i-1)*(*npat),alpha0+i*(*nclust),nu+i*(*nclust),balpha+i,nualpha+i,rho+i*(*nclust),prob+i*(*npat),nclust,npat,h_alpha0,h_nu,h_balpha,h_nualpha,h_rho,h_prob,&logscale);

    propdens1= dproposal(alpha0+i*(*nclust),nu+i*(*nclust),balpha+i,nualpha+i,rho+i*(*nclust),prob+i*(*npat),alpha0+(i-1)*(*nclust),nu+(i-1)*(*nclust),balpha+i-1,nualpha+i-1,rho+(i-1)*(*nclust),prob+(i-1)*(*npat),nclust,npat,h_alpha0,h_nu,h_balpha,h_nualpha,h_rho,h_prob,&logscale);



    if ((*method)==1) { //Metropolis-Hastings

      lambda= exp(lhood[i]-lhood[i-1]+propdens0-propdens1);

    } else {            //Simulated annealing

      temperature= 1.0/(log(1.0+i));

      lambda= exp((lhood[i]-lhood[i-1])/temperature);

    }

    u= runif();

    if (u<lambda) { //accept MH step

      if ((*trace ==1) && (maxlhood<lhood[i])) {

	maxlhood= lhood[i];

        if (*method ==0) printf("        %d   %.6f \n",i,lhood[i]);

      }

      (*acprop)= (*acprop)+1;

    } else {        //reject MH step

      lhood[i]= lhood[i-1];

      balpha[i]= balpha[i-1];

      nualpha[i]= nualpha[i-1];

      for (j=0; j< *nclust; j++) { 

	alpha0[i*(*nclust)+j]= alpha0[(i-1)*(*nclust)+j]; 

	nu[i*(*nclust)+j]= nu[(i-1)*(*nclust)+j]; 

	rho[i*(*nclust)+j]= rho[(i-1)*(*nclust)+j]; 

      }

      for (j=0; j< *npat; j++) { prob[i*(*npat)+j]= prob[(i-1)*(*npat)+j]; }

    }

    if ((*trace==1) && (*method==1) && ((i % B10)==0)) printf("  %d iterations \n",i);

  }

  (*acprop)= (*acprop)/(.0+ *B);

  if ((*trace==1) && (*method==1)) printf("Proportion of accepted proposals: %.2f \n",*acprop);



  free_ivector(sel,0,*nrow);

  free_dvector(sumd,0,*npat); free_dvector(sumci,0,*nclust); free_dvector(v,0,(*nrow)*(*npat));

  free_dvector(sumx,0,(*nrow)*ncolsumx); free_dvector(nobsx,0,ncolsumx);

  free_dvector(sumxpred,0,(*nrow)*ncolsumx);  free_dvector(nobsxpred,0,ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodx,0,(*nrow));

    free_dvector(prodxpred,0,(*nrow));

  } else {

    free_dvector(prodx,0,(*nrow)*ncolsumx);

    free_dvector(prodxpred,0,(*nrow)*ncolsumx);

  }



  free_ivector(cluslist,0,*nclust);

  free_dvector(ngroupstot,0,*nclust);



}





double dprior_hyperpar(double *alpha0, double *nu, double *balpha, double *nualpha, double *rho, double *prob, int *nclust, int *npat, double *a_alpha0, double *b_alpha0, double *a_nu, double *b_nu, double *a_balpha, double *b_balpha, double *a_nualpha, double *b_nualpha, double *p_rho, double *p_prob, int *logscale) {

  //Evaluates prior density of hyper-parameters

  int i; double ans;



  for (i=0, ans=0; i< *nclust; i++) {

    ans+= log(dgammaC(alpha0[i],*a_alpha0,*b_alpha0)) + log(dinvgammaC(nu[i],*a_nu,*b_nu));

  }

  ans+= log(dgammaC(*balpha,*a_balpha,*b_balpha)) + log(dinvgammaC(*nualpha,*a_nualpha,*b_nualpha));

  ans+= log(ddirichlet(prob,p_prob,npat));

  if (*nclust >1) ans+= log(ddirichlet(rho,p_rho,nclust));

  if (*logscale==1) return(ans); else return(exp(ans));



}



double dproposal(double *alpha0new, double *nunew, double *balphanew, double *nualphanew, double *rhonew, double *probnew, double *alpha0, double *nu, double *balpha, double *nualpha, double *rho, double *prob, int *nclust, int *npat, double *h_alpha0, double *h_nu, double *h_balpha, double *h_nualpha, double *h_rho, double *h_prob, int *logscale) {

/* Evaluates hyper-par proposal density (used for MH rejection prob)

     alpha0new ~ Gamma(h_alpha0,h_alpha0/alpha0)

     nunew ~ IGamma(h_nu,h_nu*nu)

     balphanew ~ Gamma(h_bapha,h_balpha/balpha)

     nualphanew ~ IGamma(h_nualpha,nualpha*h_nualpha)

     rho ~ Dirichlet(rho*h_rho)

     prob ~ Dirichlet(prob*h_prob)

   Returns: density of the proposal

*/

  int i; double *rhopar, *probpar, ans;



  rhopar= dvector(0,*nclust); probpar= dvector(0,*npat);

  for (i=0; i< *npat; i++) { probpar[i]= prob[i]*(*h_prob); }

  for (i=0, ans=0; i< *nclust; i++) {

    ans+= log(dgammaC(alpha0new[i],*h_alpha0,(*h_alpha0)/alpha0[i])) + log(dinvgammaC(nunew[i],*h_nu,(*h_nu)*nu[i]));

    rhopar[i]= rho[i]*(*h_rho);

  }

  ans+= log(dgammaC(*balphanew,*h_balpha,(*h_balpha)/(*balpha))) + log(dinvgammaC(*nualphanew,*h_nualpha,(*h_nualpha)*(*nualpha)));

  if (*nclust >1) ans+= ddirichlet(rhonew,rhopar,nclust);

  ans+= ddirichlet(probnew,probpar,npat);



  free_dvector(rhopar,0,*nclust); free_dvector(probpar,0,*npat);

  if (*logscale==1) return(ans); else return(exp(ans));



}





void rproposal(double *alpha0new, double *nunew, double *balphanew, double *nualphanew, double *rhonew, double *probnew, double *alpha0, double *nu, double *balpha, double *nualpha, double *rho, double *prob, int *nclust, int *npat, double *h_alpha0, double *h_nu, double *h_balpha, double *h_nualpha, double *h_rho, double *h_prob) {

/* Evaluates hyper-par proposal density (used for MH rejection prob)

     alpha0new ~ Gamma(h_alpha0,h_alpha0/alpha0)

     nunew ~ IGamma(h_nu,h_nu*nu)

     balphanew ~ Gamma(h_bapha,h_balpha/balpha)

     nualphanew ~ IGamma(h_nualpha,nualpha*h_nualpha)

     rho ~ Dirichlet(rho*h_rho)

     prob ~ Dirichlet(prob*h_prob)

   Returns: proposed hyper-par values in alpha0new, nunew, balphanew, nualphanew, rhonew, probnew

*/

  int i; double *rhopar, *probpar;



  rhopar= dvector(0,*nclust); probpar= dvector(0,*npat);

  for (i=0; i< *nclust; i++) { rhopar[i]= rho[i]*(*h_rho); }

  for (i=0; i< *npat; i++) { probpar[i]= prob[i]*(*h_prob); }



  for (i=0; i< *nclust; i++) {

    alpha0new[i]= rgammaC(*h_alpha0,(*h_alpha0)/alpha0[i]);

    nunew[i]= 1.0/rgammaC(*h_nu,(*h_nu)*nu[i]);

  }

  (*balphanew)= rgammaC(*h_balpha,(*h_balpha)/(*balpha));

  (*nualphanew)= 1.0/rgammaC(*h_nualpha,(*h_nualpha)*(*nualpha));

  if (*nclust >1) rdirichlet(rhonew,rhopar,nclust); else *rhonew= 1;

  rdirichlet(probnew,probpar,npat);



  free_dvector(rhopar,0,*nclust); free_dvector(probpar,0,*npat);



}



/*********************************************************************************************

                                   TERMINAL UTILITY ROUTINES

*********************************************************************************************/



void sampleclas_ggC(int *d, double *pgroup, double *xnew, int *nsel, int *sel, int *nrow, int *ncol, double *x, int *groups, int *K, double *Kprob, double *rho, double *prob, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox) {

/* Computes posterior prob that a new sample belongs to each group and classifies it to the group with highest prob

   Input arguments

   - xnew: vector of length nsel containing the observations for the new sample

   - nsel: number of observations to be used for the classification

   - sel: vector of length nsel indicating the row numbers that the observations in xnew correspond to in v, sumx, and nobsx.

   - nrow: number of rows of x (genes). 

   - ncol: number of columns of x (samples).

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - Kprob: vector with prior probabilities for each group

   - rho: vector with estimated prob of each cluster (length *nclust)

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - alpha0, nu: hyper-prior for lambda is IGamma(alpha0,alpha0/nu)

   - balpha, nualpha: hyper-prior for alpha is Gamma(balpha,balpha/nualpha)

   - equalcv: set to 1 to indicate CV is constant across groups, 0 otherwise

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on.

   - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

   - nobsx: number of terms in the sum. As for sumx, the first ngrouppat[0] cols correspond to pattern 0 and so on.

   - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, nobsx. For usesumx==0 they're computed from the data x.

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used

*/

/* Output arguments

   - d: group that the new sample is classified into

   - pgroup: vector of length *K with the posterior probabilities of each group

*/



  int i, j, k, m, ncolsumx, *colini, one=1, coldi, newton=2, logscale=1, usexpred=0, lengtha;

  double logl, sumf, a2, b1, b2, *s, psum, maxp, k1, k2, k3, *sumxpred, *prodxpred, *nobsxpred;



  sumxpred= dvector(0,1); prodxpred= dvector(0,1); nobsxpred= dvector(0,1); //never used: initialized to avoid warnings



  colini= ivector(0,*npat);

  colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

  ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];



  if (*usesumx == 0) { compute_sumxC(sumx,prodx,nobsx,equalcv,nsel,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one); }



  if (*equalcv == 1) {      //CONSTANT CV ACROSS GROUPS

    psum= 0;

    for (k=0; k<(*K); k++) {                                       //loop over groups

      logl=0;

      for (i=0; i<(*nsel); i++) {                                  //loop over nsel genes

	sumf= 0;

	for (j=0; j<(*npat); j++) {                                //sum over patterns

	  coldi= colini[j];

	  for (m=0; m<(*nclust); m++) {                            //sum over clusters

	    lengtha= ngrouppat[j];

	    b1= balpha[0]/nualpha[0] - prodx[sel[i]]; b2= alpha0[m]/nu[m];

	    s= sumx + sel[i]*ncolsumx;

	    k1= kcgammaC(nobsx+coldi,balpha,&b1,alpha0+m,&b2,s+coldi,&lengtha,&newton,&logscale);

	    b1= balpha[0]/nualpha[0] - prodx[sel[i]] - log(xnew[i]);

	    s[coldi+patterns[(*K)*j+k]]+= xnew[i];

	    nobsx[coldi+patterns[(*K)*j+k]]+= 1;

	    k2= kcgammaC(nobsx+coldi,balpha,&b1,alpha0+m,&b2,s+coldi,&lengtha,&newton,&logscale);

	    s[coldi+patterns[(*K)*j+k]]-= xnew[i];

	    nobsx[coldi+patterns[(*K)*j+k]]-= 1;

	    k3= pdfcond_pat_clus(sel[i],j,m,alpha0,nu,balpha,nualpha,ngrouppat,colini,ncolsumx,sumx,sumxpred,prodx,prodxpred,nobsx,nobsxpred,equalcv,usexpred);

	    sumf += exp(k2-k1-log(xnew[i])+k3) * prob[j] * rho[m];

	  }

	} //end j loop

	logl += log(sumf);

      } //end i loop

      //un-normalized post prob

      if (k==0) { pgroup[0]= logl; psum= 1; } else { pgroup[k]= exp(logl-pgroup[0])*Kprob[k]/Kprob[0]; psum+= pgroup[k]; }

    } //end k loop

  } else {               //VARYING CV ACROSS GROUPS

    lengtha = 1;

    s= dvector(0,1);

    psum= 0;

    for (k=0; k<(*K); k++) {                                       //loop over groups

      logl=0;

      for (i=0; i<(*nsel); i++) {                                  //loop over nsel genes

	sumf= 0;

	for (j=0; j<(*npat); j++) {                                //sum over patterns

	  coldi= colini[j]+patterns[(*K)*j+k];

	  for (m=0; m<(*nclust); m++) {                            //sum over clusters

	    b1= balpha[0]/nualpha[0] - prodx[sel[i]*ncolsumx+coldi]; b2= alpha0[m]/nu[m];

	    s[0]= sumx[sel[i]*ncolsumx+coldi];

	    k1= kcgammaC(nobsx+coldi,balpha,&b1,alpha0+m,&b2,s,&lengtha,&newton,&logscale);

	    b1= balpha[0]/nualpha[0] - prodx[sel[i]*ncolsumx+coldi] - log(xnew[i]);

	    s[0]= sumx[sel[i]*ncolsumx+coldi] + xnew[i];

	    a2= nobsx[coldi]+1;

	    k2= kcgammaC(&a2,balpha,&b1,alpha0+m,&b2,s,&lengtha,&newton,&logscale);

	    k3= pdfcond_pat_clus(sel[i],j,m,alpha0,nu,balpha,nualpha,ngrouppat,colini,ncolsumx,sumx,sumxpred,prodx,prodxpred,nobsx,nobsxpred,equalcv,usexpred);

	    sumf += exp(k2-k1-log(xnew[i])+k3) * prob[j] * rho[m];

	  }

	} //end j loop

	logl += log(sumf);

      } //end i loop

      //un-normalized post prob

      if (k==0) { pgroup[0]= logl; psum= 1; } else { pgroup[k]= exp(logl-pgroup[0])*Kprob[k]/Kprob[0]; psum+= pgroup[k]; }

    } //end k loop

    free_dvector(s,0,1);

  }



  maxp= pgroup[0]= 1/psum; (*d)= 0;

  for (k=1; k<(*K); k++) {                        //normalize post prob & find most probable group

    pgroup[k]= pgroup[k]/psum;

    if (pgroup[k]>maxp) { (*d)= k; maxp= pgroup[k]; }

  }



  free_dvector(sumxpred,0,1); free_dvector(prodxpred,0,1); free_dvector(nobsxpred,0,1);



}





void maxwtpfp(double *u, int *d, double *fdr, double *fnr, double *cf, int *nsel, int *sel, double *v, int *npat) {

/* Optimal terminal decisions to max cf[0]*E(TP0) - cf[1]*E(FP0) + cf[2]*E(TP1) - cf[3]*E(FP1)

   Input arguments

   cf: coefficients for the utility function

   nsel: terminal decisions are only taken for nsel genes with lowest prob of being equally expressed

   sel: indexes of the nsel genes

   v: vector with posterior probabilities. It's really a matrix with genes in rows and expression patterns in cols, entered in row order.

   npat: number of columns of v

*/

/* Ouput arguments

   u: expected utility

   d: pattern to which each gene is assigned to

   fdr: E(FDR)

   fnr: E(FNR)

*/



  int i,j,jmax,nde=0,nee=0;

  double v0=0,vmax=0,fd=0,fn=0;



  (*u)= 0;

  for (i=0;i<(*nsel);i++) {

    v0= v[sel[i]*(*npat)];                                                        //find non-null pattern with highest probability

    vmax= v[sel[i]*(*npat)+1]; jmax=1;

    for (j=2;j<(*npat);j++) {                                                    

      if (v[sel[i]*(*npat)+j]>vmax) { jmax=j; vmax= v[sel[i]*(*npat)+j]; }

    }

    if (vmax > v0*(cf[0]+cf[1])/(cf[2]+cf[3])) {                                  //determine if it's optimal to assign gene to pattern with highest prob

      d[sel[i]]= jmax;

      (*u) += cf[2]*vmax - cf[1]*v0 - cf[3]*(1-v0-vmax);

      fd+= v0; nde++;

    } else {

      d[sel[i]]=0;

      (*u) += cf[0]*v0 - cf[3]*(1-v0);

      fn+= 1-v0; nee++;

    }

  }                                                                               //end i for



  if (nde>0) { (*fdr)= fd/nde; } else { (*fdr)= 0; }                              //find E(FDR), E(FNR)

  if (nde<(*nsel)) { (*fnr) = fn/nee; } else { (*fnr)= 0; }



}





void minfnrstfdr(double *u, double *threshold, int *d, double *fdr, double *fnr, double *power, int *nsel, int *sel, double *v, int *npat, double *fdrmax) {

/* Optimal terminal decisions to min E(FNR) s.t. E(fdr)<fdrmax. */

/* Input arguments

    nsel: terminal decisions are only taken for nsel genes with lowest prob of being equally expressed

    sel: indexes of the nsel genes

    v: vector with posterior probabilities. It's really a matrix with genes in rows and expression patterns in cols, entered in row order.

    npat: number of columns of v

    fdrmax: upper bound for E(FDR)

*/

/* Output arguments

    u: expected number of True Positives

    threshold: optimal threshold for posterior prob of pattern 0 (genes with prob < threshold are declared DE)

    d: pattern to which each gene is assigned to

    fdr: E(FDR)

    fnr: E(FNR)

    power: power as estimated by E(TP)/E(Positives)

*/



int i, j, jopt, nde, nee;

double *v0, vopt, fn, fd, ntrue;



v0 = dvector(0,(*nsel));                                                         //put prob of pattern 0 into v0

for (i=0;i<(*nsel);i++) { v0[i]= v[sel[i]*(*npat)]; }

dvecsort(v0,(*nsel));                                                            //sort v0 in ascending order



fd= fn= ntrue= 0; nde= nee= 0; (*threshold)= 0;



for (i=0; i<(*nsel); i++) {                                                      //find the optimal threshold

  if ((fd+v0[i])/(nde+1) < (*fdrmax)) {

    (*threshold)= v0[i];

    fd += v0[i]; nde++;

  } else {

    fn += 1-v0[i]; nee++;

  }

  ntrue+= 1-v0[i];

}                                                                               //end i for



if (nde>0) { (*fdr)= fd/nde; } else { (*fdr)= 0; }                              //find E(FDR), E(FNR)

if (nde<(*nsel)) { (*fnr) = fn/nee; } else { (*fnr)= 0; }



(*u) = 0;

for (i=0; i<(*nsel); i++) {                                                     //assign genes to expression patterns

  if (v[sel[i]*(*npat)] <= (*threshold)) {

    jopt=1; vopt= v[sel[i]*(*npat)+1];

    for (j=1; j<(*npat); j++) { if (v[sel[i]*(*npat)+j]>v[sel[i]*(*npat)+jopt]) { jopt= j; vopt= v[sel[i]*(*npat)+j]; } }

    d[sel[i]]= jopt;

    (*u) += vopt;

  } else {

    d[sel[i]]= 0;

  }

}                                                                               //end i for



(*power)= (*u)/ntrue;



free_dvector(v0,0,(*nsel));



}



void geneclus(int *d, double *ppat, int *nsel, int *sel, double *v, int *npat) {

/* Classifies each gene into the highest posterior probability pattern */

/* Input arguments

    nsel: terminal decisions are only taken for nsel genes with lowest prob of being equally expressed

    sel: indexes of the nsel genes

    v: vector with posterior probabilities. It's really a matrix with genes in rows and expression patterns in cols, entered in row order.

    npat: number of columns of v

*/

/* Output arguments

    d: pattern to which each gene is assigned to

    ppat: posterior probability of the assigned pattern

*/



int i, j, jopt;



  for (i=0; i<(*nsel); i++) {                          

    for (j=1, jopt=0; j<(*npat); j++) { if (v[sel[i]*(*npat)+j]>v[sel[i]*(*npat)+jopt]) jopt= j; }

    d[sel[i]]= jopt;

    ppat[sel[i]]= v[sel[i]*(*npat)+jopt];

  }                                                    



}





void sel_mostDEgenes(int *nsel, int *sel, int *genelimit, double *v0thre, int *nrow, int *npat, double *v) {

/* Selects the nsel genes with the lowest probability of being equally expressed */

/* Input

   - genelimit: only the genelimit genes with lowest v0 (prob of being equally expressed across all groups) will be used to classify samples

   - v0thre: only the genes with v0thre>=v0 (prob of being equally expressed across all groups) are used to classify samples. If no genes make the threshold, the classification is based on the gene with the smallest v0.

   - nrow: number of genes

   - npat: number of expression patterns

   - v: matrix (nrow rows, npat cols) with posterior probabilities for each gene and expression pattern

*/

/* Ouput

   - nsel: number of genes selected

   - sel: indexes of genes that were selected (only positions 0:nsel-1 are updated).

 */

  int i, selmin;

  double v0min, *v0sel;



  v0sel= dvector(0,*nrow);



  selmin= 0; v0min= v[0]; *nsel=0;

  for (i=0; i<(*nrow); i++) {                                                             //select genes with v0<=v0thre and find min(v0)

    if (v[i*(*npat)] <= (*v0thre)) { v0sel[*nsel]= v[i*(*npat)]; sel[*nsel]= i; (*nsel)++; }

    if (v[i*(*npat)] < v0min) { v0min= v[i*(*npat)]; selmin= i; }

  }

  if (*nsel==0) {                                                                          //if no genes selected, select gene with smallest v0

    v0sel[0]= v0min; sel[0]= selmin; (*nsel)= 1;

  } else if ((*nsel)>(*genelimit)) {                                                       //if too many selected, keep genelimit genes with smallest v0

    dindexsort(v0sel,sel,0,*nsel-1,1);

    *nsel= *genelimit;

  }



  free_dvector(v0sel,0,*nrow);



}





void uobsgeneC(double *u, double *fdr, double *fnr, int *util, int *nsel, int *sel, int *d, int *trued, double *cf) {

/* Realized terminal utility, FDR & FNR for simulated data (i.e. assuming true expression pattern is known)

   Input arguments

   nsel: terminal decisions are only taken for nsel genes with lowest prob of being equally expressed

   sel: indexes of the nsel genes

   util: util==1 means utility in minfnrstfdr, util==2 maxwtpfp

   d: pattern to which each gene is assigned to (typically, as returned by utgene_parC)

   trued: true differential expression pattern

   cf: coefficients for the utility function. Ignored if util==1.

   nrow: number of genes

*/

/* Output arguments

   u: realized utility (for util==1 returns number of true positives)

   fdr: realized FDR

   fnr: realized FNR

*/



  int i;

  double fd=0, fn=0, nde=0, nee=0;



  (*u)= 0;

  for (i=0; i<(*nsel); i++) {

    if (d[sel[i]]!=0) { 

      if (trued[sel[i]]==0) fd++; 

      nde++;

    } else {

      if (trued[sel[i]]!=0) fn++; 

      nee++;

    }

    if ((*util)==1) {

      (*u)+= ((d[sel[i]]==trued[sel[i]]) && (trued[sel[i]]!=0));

    } else if ((*util)==2) {

      (*u)+= cf[0]*(d[sel[i]]==trued[sel[i]] && trued[sel[i]]==0) - cf[1]*((d[sel[i]]!=trued[sel[i]]) && trued[sel[i]]==0) + cf[2]*(d[sel[i]]==trued[sel[i]] && trued[sel[i]]!=0) - cf[3]*((d[sel[i]]!=trued[sel[i]]) && trued[sel[i]]!=0);

    }

  }                                                                               //end i for



  if (nde>0) { (*fdr)= fd/nde; } else { (*fdr)= 0; }                              //find FDR, FNR

  if (nde<(*nsel)) { (*fnr) = fn/nee; } else { (*fnr)= 0; }



}



void expected_fp(double *efp, double *fdrseq, int *nthre, int *B, int *niter, double *z, double *m, double *s, int *index, int *znclust, int *zclustsize, int *nrow, int *ncol, double *x, int *groups, int *K, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, int *cluslist, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, double *sumxpred, double *prodxpred, double *nobsxpred, int *gapprox) {

/* Estimate expected nb of false positives in gene DE analysis by permutation methods for

   a range of fdrseqs on the posterior prob of EE */

/* Input arguments

   - fdrseq: sequence of thresholds for target FDR

   - nthre: length of threshold

   - B: number of permutations

   - niter: number of iterations to re-estimate prob of each expr pattern

   - z: matrix with z-scores to use for simulation under the null

   - m: mean of each row in x

   - s: sd of each row in x

   - index: indicates ordering of genes such that genes from 1st cluster come before 2nd cluster etc.

   - znclust: number of clusters for the z scores

   - zclustsize: number of genes (row of x) in each cluster (vector of length znclust)

   - nrow: number of genes (number of rows in x)

   - other arguments passed to pp_ggC. Careful: rows of x must be ordered 

   Output arguments

   - efp: mean number of false positives for each fdrseq, averaged over B permutations

*/

  int *groups0, i, j, k, *sel, usesumx=0, *nde, *nde0, B10, u;

  double *v, lhood, *prob0, *probnew, *xboot, *threshold, *pboot;

  //FILE *ifile; //debug



  if (*B < 10) (*B)= 10; //ensure a minimum amount of permutations

  nde= ivector(0,*nthre); nde0= ivector(0,*nthre); threshold= dvector(0,*nthre);

  sel= ivector(0,*nrow);

  for (i=0; i<(*nrow); i++) { sel[i]= i; }

  groups0= ivector(0,*ncol);

  for (i=0; i<(*ncol); i++) { groups0[i]= groups[i]; }

  probnew=dvector(0,*npat); prob0= dvector(0,*npat);

  for (i=0; i<(*npat); i++) { prob0[i]= prob[i]; }

  v= dvector(0, (*nrow)*(*npat));

  dvecsort(fdrseq,*nthre);  //put sequence of fdrseq in ascending order

  for (i=0; i<(*nthre); i++) { efp[i]= 0; }

  pp_ggC(v,&lhood,nrow,sel,ncol,x,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob0,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&usesumx,gapprox); //post prob

  countde(nde,threshold,nthre,fdrseq,nrow,v,npat);   //count nb of genes declared DE

  xboot= dvector(0,(*nrow)*(*ncol)); 

  pboot= dvector(0,*ncol);

  for (i=0; i<(*ncol); i++) { pboot[i]= 1.0/(*ncol); }



  for (i=0, B10=(*B/10); i<(*B); i++) {   //loop over permutations

    if (*znclust > 0) {

      bootnull(xboot,nrow,ncol,z,m,s,index,znclust,zclustsize);  //gene cluster level bootstrap

    } else {

      for (k=0; k<(*nrow); k++) {

	for (j=0; j<(*ncol); j++) {

          u= runifdisc(0,*ncol - 1);

	  xboot[k*(*ncol)+j]= x[k*(*ncol)+u];

	}

      }

    }



    //for (j=0,err=1; (j<(*niter)) && (err>eps); j++) {     //re-estimate prob of each expr pattern

    //pp_ggC(v,&lhood,nrow,sel,ncol,xboot,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob0,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&usesumx,gapprox); //post prob

    //colMeans(probnew,v,*nrow,*npat);

    //for (err=0,k=0; k<(*npat); k++) { err += fabs(probnew[k]-prob0[k]); prob0[k]= probnew[k]; }

    //}

    //countde(nde0,nthre,fdrseq,nrow,v,npat);   //nb of DE genes re-applying optimal rule

    pp_ggC(v,&lhood,nrow,sel,ncol,xboot,groups,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob0,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&usesumx,gapprox); //pp wo re-estimation

    countde_threshold(nde0,threshold,nthre,nrow,v,npat); //nb of DE genes using same rejection region



    for (j=0; j<(*nthre); j++) { efp[j] += nde0[j]; }  //add to nb of false positives

    //if (((i+1) % B10)==0) printf("  %d done \n",i+1);

  }

  for (i=0; i<(*nthre); i++) {        //average

    if (nde[i]>0) { efp[i] = prob[0]*efp[i]/((*B+.0)*(nde[i]+.0)); } else { efp[i]= 0; }

  } 



  free_ivector(nde,0,*nthre); free_ivector(nde0,0,*nthre); free_dvector(threshold,0,*nthre);

  free_ivector(sel,0,*nrow);

  free_ivector(groups0,0,*ncol);

  free_dvector(probnew,0,*npat); free_dvector(prob0,0,*npat);

  free_dvector(v,0, (*nrow)*(*npat));

  free_dvector(xboot,0,(*nrow)*(*ncol)); free_dvector(pboot,0,*ncol);



}



void bootnull(double *xboot, int *nrow, int *ncol, double *z, double *m, double *s, int *index, int *znclust, int *zclustsize) {

//Obtain bootstrap sample under the null, by sampling z-scores from genes in the same cluster

/* Input arguments

   - nrow, ncol: number of rows and columns of xboot

   - z: matrix with z-scores to use for simulation under the null

   - m: mean of each row in x

   - s: sd of each row in x

   - index: indicates ordering of genes such that genes from 1st cluster come before 2nd cluster etc.

   - znclust: number of clusters for the z scores

   - zclustsize: number of genes (row of x) in each cluster (vector of length znclust)

   Output arguments

   - xboot: bootstrap sample taken by selecting z-scores from genes in the same cluster, then scaling by s and shifting by m

*/

  int i, j, k, *posmin, *posmax, u, ucol, nrep;



  posmin= ivector(0,*znclust); posmax= ivector(0,*znclust);

  posmin[0]= 0; posmax[0]= zclustsize[0]-1;

  for (i=1; i<(*znclust); i++) { posmin[i]= posmax[i-1]+1; posmax[i]= posmax[i-1]+zclustsize[i]; }



  for (i=0, k=0; i<(*nrow); i++) {

    if (i==(posmax[k]+1)) k++;

    for (j=0; j<(*ncol); j++) {

      u= runifdisc(posmin[k],posmax[k]); ucol= runifdisc(0,*ncol-1);

      xboot[index[i]*(*ncol)+j]= z[index[u]*(*ncol)+ucol]*s[index[i]] + m[index[i]];

      nrep= 0;

      while(xboot[index[i]*(*ncol)+j]<0) {  //ensure positivity

        u= runifdisc(posmin[k],posmax[k]); ucol= runifdisc(0,*ncol-1);

        xboot[index[i]*(*ncol)+j]= z[index[u]*(*ncol)+ucol]*s[index[i]] + m[index[i]];

	if (nrep==5) { xboot[index[i]*(*ncol)+j]= m[index[i]]; } else { nrep++; }

      }

    }

  }

  free_ivector(posmin,0,*znclust); free_ivector(posmax,0,*znclust);



}





void countde_threshold(int *nde, double *threshold, int *nthre, int *nrow, double *v, int *npat) {

//Count nb of genes with post prob of EE below a certain threshold

/* Input arguments

   - threshold: sequence of thresholds

   - nthre: length of threshold

   - nrow: number of genes (rows in v)

   - v: matrix with posterior prob of each expression pattern (nrow rows, npat cols)

   - npat: number of expression patterns (cols in v)

   Output arguments

   - nde: number of genes declared DE for each threshold

*/



  int i, j;

  double *v0;



  v0= dvector(0,*nrow);   //copy prob of EE into v0 and put in increasing order

  for (i=0; i<(*nrow); i++) { v0[i]= v[i*(*npat)]; }

  dvecsort(v0,*nrow);



  for (i=0; i<(*nthre); i++) { nde[i]= 0; }

  for (i=0, j=0; (i<(*nrow)) && (j<(*nthre)); i++) {

    while ((v0[i]>threshold[j]) && (j<(*nthre))) j++; //FDR at which gene rejected

    if ((v0[i]<=threshold[j]) && (j<(*nthre))) {   //if gene indeed rejected

      nde[j] +=1;

    }

  }

  for (i=1; i<(*nthre); i++) { nde[i] += nde[i-1]; }



  free_dvector(v0,0,*nrow);



}





void countde(int *nde, double *threshold, int *nthre, double *fdrseq, int *nrow, double *v, int *npat) {

//Count nb of genes declared DE by Mueller's optimal procedure and the corresponding threshold on the posterior probability of EE, for a sequence of bounds on the E(FDR)

/* Input arguments

   - nthre: length of the sequence of bounds

   - fdrseq: sequence of bounds on the E(FDR)

   - nrow: number of genes (rows in v)

   - v: matrix with posterior prob of each expression pattern (nrow rows, npat cols)

   - npat: number of expression patterns (cols in v)

   Output arguments

   - nde: number of genes declared DE for each bound in fdrseq

   - threshold: threshold for the post prob corresponding to each target E(FDR) in fdrseq

*/



  int i, j, sumnde;

  double *v0, fd;



  v0= dvector(0,*nrow);   //copy prob of EE into v0 and put in increasing order

  for (i=0; i<(*nrow); i++) { v0[i]= v[i*(*npat)]; }

  dvecsort(v0,*nrow);



  for (i=0; i<(*nthre); i++) { nde[i]= 0; threshold[i]= 0; }

  fd= 0; sumnde= 0;

  for (i=0, j=0; (i<(*nrow)) && (j<(*nthre)); i++) {

    while (((fd+v0[i])/(sumnde+1.0) >= fdrseq[j]) && (j<(*nthre))) j++; //FDR at which gene rejected

    if (((fd+v0[i])/(sumnde+1.0) < fdrseq[j]) && (j<(*nthre))) {   //if gene indeed rejected

      fd += v0[i]; nde[j] += 1; sumnde++; threshold[j]= v0[i];

    }

  }

  for (i=1; i<(*nthre); i++) { nde[i] += nde[i-1]; }



  free_dvector(v0,0,*nrow);



}



void utgene_parC(double *u, int *d, double *fdr, double *fnr, double *power, double *threshold, int *util, double *cf, int *nsel, int *sel, double *v, int *npat, double *fdrmax) {

/* Optimal terminal decisions and expected utility for several gene differential expression analysis utilities

   Input arguments

   util: util==1 means calling minfnrstfdr, util==2 calls maxwtpfp.

   cf: coefficients for the utility function. Ignored if util==1.

   nsel: terminal decisions are only taken for nsel genes with lowest prob of being equally expressed

   sel: indexes of the nsel genes

   v: vector with posterior probabilities. It's really a matrix with genes in rows and expression patterns in cols, entered in row order.

   npat: number of columns of v

   fdrmax: upper bound for E(FDR). Ignored if util!=1.

*/

/* Ouput arguments

   u: expected utility (for util==1 returns expected number of True Positives)

   d: pattern to which each gene is assigned to

   fdr: E(FDR)

   fnr: E(FNR)

   power: power as measured by E(TP)/E(positives)

   threshold: If util==1, optimal threshold for posterior prob of pattern 0 (genes with prob < threshold are declared DE). Ignored if util!=1.

*/



  if ((*util)==1) { 

    minfnrstfdr(u, threshold, d, fdr, fnr, power, nsel, sel, v, npat, fdrmax);

  } else if ((*util)==2) {

    maxwtpfp(u, d, fdr, fnr, cf, nsel, sel, v, npat);

  }



}





void utsample_ggC(double *ccall, double *seccall, double *ccgroup, int *ngroup, int *B, double *preceps, int *genelimit, double *v0thre, int *nsel, int *sel, int *usesel, int *nrow, int *ncol, double *x, int *groups, double *v, int *K, double *Kprob, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, int *ncolsumx, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox) {

/* Estimates expected terminal utility for sample classification. Currenly utility defined as correct classification

   Input:

   - B: maximum number of MC simulations

   - preceps: procedure stops early if .5*length of CI < preceps, but never before a minimum of 50 simulations

   - genelimit: only the genelimit genes with lowest v0 (prob of being equally expressed across all groups) will be used to classify samples

   - v0thre: only the genes with v0thre>=v0 (prob of being equally expressed across all groups) are used to classify samples. If no genes make the threshold, the classification is based on the gene with the smallest v0.

   - nsel: if usesel==1 expression values are simulated for nsel genes only and genelimit, v0thre are ignored. If usesel==0 nsel is ignored.

   - sel: if usesel==1 vector of length nrow containing the keys of the genes to simulate for. If usesel==0 it's determined based on v0thre and genelimit and returned

   - usesel: if usesel==1 expression values are simulated for genes indicated in the sel vector. If usesel==0 genelimit and v0thre are used to define what genes to simulate for

   - nrow: number of rows of x (genes). 

   - ncol: number of columns of x (samples).

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to

   - v: matrix with posterior probabilities (genes in rows, expression patterns in cols)

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - Kprob: vector with prior probabilities for each group

   - alpha0, nu: hyper-prior for lambda is IGamma(alpha0,alpha0/nu)

   - balpha, nualpha: hyper-prior for alpha is Gamma(balpha,balpha/nualpha)

   - nclust: number of clusters

   - rho: vector with estimated probability of each cluster

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - ncolsumx: number of columns of sumx

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on.

   - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

   - nobsx: number of terms in the sum. As for sumx, the first ngrouppat[0] cols correspond to pattern 0 and so on.

   - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, nobsx. For usesumx==0 they're computed from the data x.

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used

   Output:

   - ccall: estimated overall probability of correctly classifying a sample

   - seccall: standard error of the ccall estimate as estimated by SD/sum(groupsim)

   - ccgroup: number of correctly classified samples from each group

   - ngroup: number of samples from each group

   Details: ccall is estimated by using posterior probabilities of correct classification, whereas ccgroup simply counts the proportion of samples of each group that are classified correctly. ccall is therefore a better estimate than what we obtain by simply averaging the elements in ccgroup.

*/



  int i, curgroup, d, one=1, *dnew;

  double *xnew, *pgroup, *anew, *lnew, cv;



  pgroup= dvector(0,*K);



  if (*usesel == 0) { sel_mostDEgenes(nsel,sel,genelimit,v0thre,nrow,npat,v); }  //select most DE genes

  if (*usesumx == 0) { compute_sumxC(sumx,prodx,nobsx,equalcv,nsel,sel,ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one); }



  xnew= dvector(0,*nsel); dnew= ivector(0,*nsel); anew= dvector(0,*nsel); lnew= dvector(0,*nsel);

  (*ccall)= (*seccall)= 0;  

  for (i=0; i<(*K); i++) { ccgroup[i]= 0; ngroup[i]= 0; }



  cv= 0;

  for (i=0; i<(*B) && (cv>(*preceps) || i<50); i++) {                      //loop over simulations

    curgroup= rdisc(Kprob,*K);

    simnewsamples_ggC(xnew,dnew,anew,lnew,&one,&curgroup,K,nsel,sel,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,v,npat,patterns,ngrouppat,sumx,prodx,nobsx,gapprox); //get new sample

    sampleclas_ggC(&d,pgroup,xnew,nsel,sel,nrow,ncol,x,groups,K,Kprob,rho,prob,alpha0,nu,balpha,nualpha,equalcv,nclust,npat,patterns,ngrouppat,sumx,prodx,nobsx,&one,gapprox); //classify it

    (*ccall) += pgroup[d]; (*seccall) += pgroup[d]*pgroup[d];

    if (d==curgroup) { ccgroup[curgroup] += 1; }

    ngroup[curgroup] += 1;

    if (i>=50) { cv= sqrt((*seccall)/(i+1) - (*ccall)*(*ccall)/((i+1)*(i+1)))/sqrt(i+i); }

  } //end i for



  (*ccall)= (*ccall)/i;                                          //mean prob of correct classification

  (*seccall)= sqrt((*seccall)/i - (*ccall)*(*ccall))/sqrt(i);    //standard error



  free_dvector(xnew,0,*nsel); free_ivector(dnew,0,*nsel); free_dvector(anew,0,*nsel); free_dvector(lnew,0,*nsel);

  free_dvector(pgroup,0,*K);



}





/*********************************************************************************************

                              PREDICTIVE TERMINAL UTILITY ROUTINES

*********************************************************************************************/



void simnewPP_ggC(double *pp, int *deltat, int *B, double *preceps, int *util, double *cf, int *genelimit, double *v0thre, int *nsel, int *sel, int *usesel, int *nrow, int *ncol, double *x, int *groups, int *K, double *v, double *alpha0, double *nu, double *balpha, double *nualpha,  int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *fdrmax, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox) {

/* Simulate posterior probabilities from the posterior predictive distribution */

/* Input:

   - deltat: expected predictive terminal utility will be evaluated after taking deltat observations from the predictive distribution

   - B: number of MC simulations used to estimate the expected terminal utility. Procedure stops early if coefficient of variation < .01 (minimum 100 sims)

   - util: util==1 means calling minfnrstfdr, util==2 calls maxwtpfp, util==3 calls maxamountde

   - cf: coefficients for the utility function. Ignored if util==1.

   - genelimit: only the genelimit genes with lowest v0 (prob of being equally expressed across all groups) will be used to classify samples

   - v0thre: only the genes with v0thre>=v0 (prob of being equally expressed across all groups) are used to classify samples. If no genes make the threshold, the classification is based on the gene with the smallest v0.

   - nsel: if usesel==1 expression values are simulated for nsel genes only and genelimit, v0thre are ignored. If usesel==0 nsel is ignored.

   - sel: if usesel==1 vector of length nrow containing the keys of the genes to simulate for. If usesel==0 it's determined based on v0thre and genelimit and returned

   - usesel: if usesel==1 expression values are simulated for genes indicated in the sel vector. If usesel==0 genelimit and v0thre are used to define what genes to simulate for

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - v: vector with posterior probabilities. It's really a matrix with genes in rows and expression patterns in cols, entered in row order.

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - balpha: estimate for b parameter in hyper-prior for alpha

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - fdrmax: upper bound for E(FDR). Ignored if util!=1.

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on.

   - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

   - nobsx: number of terms in the sum. As for sumx, the first ngrouppat[0] cols correspond to pattern 0 and so on.

   - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, nobsx. For usesumx==0 they're computed from the data x.

   - gapprox: for gapprox==1 the normalization constant for the conjugate Gamma shape distributions is computed approximately



 Output:

   - pp: matrix with posterior probabilities of equal expression

*/





  int i, *dpred, *d, *groupspred, ncolpred, ncolsumx, zero=0, one=1, *cluslist;

  double *xpred, *apred, *lpred, *vpred, *sumxpred, *prodxpred, *nobsxpred, lhood;



  ncolpred= (*deltat)*(*K);

  xpred= dvector(0,(*nrow)*ncolpred);

  dpred= ivector(0,(*nrow));

  d= ivector(0,(*nrow));

  if (*equalcv == 1) apred= dvector(0,*nrow); else apred= dvector(0,(*nrow)*(*K));

  lpred= dvector(0,(*nrow)*(*K));

  vpred= dvector(0,(*nrow)*(*npat));

  cluslist= ivector(0,*nclust);



  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;

  groupspred= ivector(0,ncolpred);                       //indicate group each new obs will belong to

  for (i=0; i<ncolpred; i++) { groupspred[i]= floor(i / (*deltat)); }



  if (*usesel == 0) { sel_mostDEgenes(nsel,sel,genelimit,v0thre,nrow,npat,v); }  //select most DE genes



  ncolsumx= 0; for (i=0; i<(*npat); i++) { ncolsumx += ngrouppat[i]; }

  if ((*usesumx)==0) { compute_sumxC(sumx,prodx,nobsx,equalcv,nsel,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one); }

  sumxpred= dvector(0,(*nrow)*ncolsumx); nobsxpred= dvector(0,ncolsumx);

  if (*equalcv == 1) {

    prodxpred= dvector(0,(*nrow));

  } else {

    prodxpred= dvector(0,(*nrow)*ncolsumx);

  }



  for (i=0; i<(*B); i++) {                 //loop over simulations

    simpred_ggC(xpred,dpred,apred,lpred,&zero,deltat,groups,K,nsel,sel,nrow,ncol,x,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,v,npat,patterns,ngrouppat,sumx,prodx,nobsx,&one,gapprox);

    compute_sumxC(sumxpred,prodxpred,nobsxpred,equalcv,nsel,sel,&ncolsumx,&ncolpred,xpred,groupspred,K,npat,patterns,ngrouppat,&one); //update suff statistics

    pp_ggC(vpred,&lhood,nsel,sel,&ncolpred,xpred,groupspred,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox); //post prob of patterns

  }                                                                    //end loop over simulations



  free_ivector(cluslist,0,*nclust);

  free_dvector(xpred,0,(*nrow)*ncolpred);

  free_ivector(dpred,0,(*nrow));

  free_ivector(d,0,(*nrow));

  if (*equalcv == 1) free_dvector(apred,0,*nrow); else free_dvector(apred,0,(*nrow)*(*K));

  free_dvector(lpred,0,(*nrow)*(*K));

  free_dvector(vpred,0,(*nrow)*(*npat));

  free_ivector(groupspred,0,ncolpred);

  free_dvector(sumxpred,0,(*nrow)*ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodxpred,0,(*nrow));

  } else {

    free_dvector(prodxpred,0,(*nrow)*ncolsumx);

  }

  free_dvector(nobsxpred,0,ncolsumx);



}





void utgene_predC(double *m, double *s, int *deltat, int *B, double *preceps, int *util, double *cf, int *genelimit, double *v0thre, int *nsel, int *sel, int *usesel, int *nrow, int *ncol, double *x, int *groups, int *K, double *v, double *alpha0, double *nu, double *balpha, double *nualpha,  int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *fdrmax, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox) {

/* Estimates expected predictive terminal utility for gene differential expression analysis after drawing deltat more observations per group */

/* Input:

   - deltat: expected predictive terminal utility will be evaluated after taking deltat observations from the predictive distribution

   - B: number of MC simulations used to estimate the expected terminal utility. Procedure stops early if coefficient of variation < .01 (minimum 100 sims)

   - util: util==1 means calling minfnrstfdr, util==2 calls maxwtpfp, util==3 calls maxamountde

   - cf: coefficients for the utility function. Ignored if util==1.

   - genelimit: only the genelimit genes with lowest v0 (prob of being equally expressed across all groups) will be used to classify samples

   - v0thre: only the genes with v0thre>=v0 (prob of being equally expressed across all groups) are used to classify samples. If no genes make the threshold, the classification is based on the gene with the smallest v0.

   - nsel: if usesel==1 expression values are simulated for nsel genes only and genelimit, v0thre are ignored. If usesel==0 nsel is ignored.

   - sel: if usesel==1 vector of length nrow containing the keys of the genes to simulate for. If usesel==0 it's determined based on v0thre and genelimit and returned

   - usesel: if usesel==1 expression values are simulated for genes indicated in the sel vector. If usesel==0 genelimit and v0thre are used to define what genes to simulate for

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - v: vector with posterior probabilities. It's really a matrix with genes in rows and expression patterns in cols, entered in row order.

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - balpha: estimate for b parameter in hyper-prior for alpha

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - fdrmax: upper bound for E(FDR). Ignored if util!=1.

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on.

   - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

   - nobsx: number of terms in the sum. As for sumx, the first ngrouppat[0] cols correspond to pattern 0 and so on.

   - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, nobsx. For usesumx==0 they're computed from the data x.

   - gapprox: for gapprox==1 the normalization constant for the conjugate Gamma shape distributions is computed approximately



 Output:

   - m: estimated expected terminal utility (sample mean of B simulations)

   - s: estimate SE, i.e. SD of the simulations / sqrt(number of simulations)

*/





  int i, *dpred, *d, *groupspred, ncolpred, ncolsumx, zero=0, one=1, *cluslist;

  double *xpred, *apred, *lpred, *vpred, uobs, fdrobs, fnrobs, powobs, threshold, *sumxpred, *prodxpred, *nobsxpred, cv, lhood;



  ncolpred= (*deltat)*(*K);

  xpred= dvector(0,(*nrow)*ncolpred);

  dpred= ivector(0,(*nrow));

  d= ivector(0,(*nrow));

  if (*equalcv == 1) apred= dvector(0,*nrow); else apred= dvector(0,(*nrow)*(*K));

  lpred= dvector(0,(*nrow)*(*K));

  vpred= dvector(0,(*nrow)*(*npat));

  cluslist= ivector(0,*nclust);



  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;

  groupspred= ivector(0,ncolpred);                       //indicate group each new obs will belong to

  for (i=0; i<ncolpred; i++) { groupspred[i]= floor(i / (*deltat)); }



  if (*usesel == 0) { sel_mostDEgenes(nsel,sel,genelimit,v0thre,nrow,npat,v); }  //select most DE genes



  ncolsumx= 0; for (i=0; i<(*npat); i++) { ncolsumx += ngrouppat[i]; }

  if ((*usesumx)==0) { compute_sumxC(sumx,prodx,nobsx,equalcv,nsel,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one); }

  sumxpred= dvector(0,(*nrow)*ncolsumx); nobsxpred= dvector(0,ncolsumx);

  if (*equalcv == 1) {

    prodxpred= dvector(0,(*nrow));

  } else {

    prodxpred= dvector(0,(*nrow)*ncolsumx);

  }



  (*m)=0; (*s)=0; cv= 1;

  for (i=0; i<(*B) && (cv>(*preceps) || i<100); i++) {                 //loop over simulations

    simpred_ggC(xpred,dpred,apred,lpred,&zero,deltat,groups,K,nsel,sel,nrow,ncol,x,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,v,npat,patterns,ngrouppat,sumx,prodx,nobsx,&one,gapprox);

    compute_sumxC(sumxpred,prodxpred,nobsxpred,equalcv,nsel,sel,&ncolsumx,&ncolpred,xpred,groupspred,K,npat,patterns,ngrouppat,&one); //update suff statistics

    pp_ggC(vpred,&lhood,nsel,sel,&ncolpred,xpred,groupspred,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox); //post prob of patterns

    utgene_parC(&uobs, d, &fdrobs, &fnrobs, &powobs, &threshold, util, cf, nsel, sel, vpred, npat, fdrmax);  //optimal term decisions & exp utility

    (*m) += uobs; (*s) += uobs * uobs;

    if ((*m)>0) { cv = sqrt((*s) / (i+1) - (*m)/(i+1) * (*m)/(i+1)) * sqrt(i+1) / (*m); } else { cv = 1; }

  }                                                                    //end loop over simulations

  (*m)= (*m) / i;

  (*s)= sqrt((*s) / i - (*m)*(*m)) / sqrt(i);



  free_ivector(cluslist,0,*nclust);

  free_dvector(xpred,0,(*nrow)*ncolpred);

  free_ivector(dpred,0,(*nrow));

  free_ivector(d,0,(*nrow));

  if (*equalcv == 1) free_dvector(apred,0,*nrow); else free_dvector(apred,0,(*nrow)*(*K));

  free_dvector(lpred,0,(*nrow)*(*K));

  free_dvector(vpred,0,(*nrow)*(*npat));

  free_ivector(groupspred,0,ncolpred);

  free_dvector(sumxpred,0,(*nrow)*ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodxpred,0,(*nrow));

  } else {

    free_dvector(prodxpred,0,(*nrow)*ncolsumx);

  }

  free_dvector(nobsxpred,0,ncolsumx);



}



void utsample_predC(double *ccall, double *seccall, double *ccgroup, int *ngroup, int *deltat, int *B, double *preceps, int *genelimit, double *v0thre, int *nsel, int *sel, int *usesel, int *nrow, int *ncol, double *x, int *groups, int *K, double *Kprob, double *v, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox) {

/* Estimates expected predictive terminal utility for sample classification after drawing deltat more observations per group */

/* Input:

   - deltat: expected predictive terminal utility will be evaluated after taking deltat observations from the predictive distribution

   - B: number of MC simulations used to estimate the expected terminal utility. 

   - preceps: procedure stops early if .5*length of CI < preceps, but never before a minimum of 100 simulations

   - genelimit: only the genelimit genes with lowest v0 (prob of being equally expressed across all groups) will be used to classify samples

   - v0thre: only the genes with v0thre>=v0 (prob of being equally expressed across all groups) are used to classify samples. If no genes make the threshold, the classification is based on the gene with the smallest v0.

   - nsel: if usesel==1 expression values are simulated for nsel genes only and genelimit, v0thre are ignored. If usesel==0 nsel is ignored.

   - sel: if usesel==1 vector of length nrow containing the keys of the genes to simulate for. If usesel==0 it's determined based on v0thre and genelimit and returned

   - usesel: if usesel==1 expression values are simulated for genes indicated in the sel vector. If usesel==0 genelimit and v0thre are used to define what genes to simulate for

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: vector with observations used to fit the model. It's really a matrix with genes in rows and samples in cols, entered in row order.

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - Kprob: vector with prior probabilities for each group

   - v: vector with posterior probabilities. It's really a matrix with genes in rows and expression patterns in cols, entered in row order.

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - balpha: estimate for b parameter in hyper-prior for alpha

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on.

   - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

   - nobsx: number of terms in the sum. As for sumx, the first ngrouppat[0] cols correspond to pattern 0 and so on.

   - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, nobsx. For usesumx==0 they're computed from the data x.

   - gapprox: for gapprox==1 the normalization constant for the conjugate Gamma shape distributions is computed approximately



 Output:

   - ccall: estimated overall probability of correctly classifying a sample

   - seccall: standard error of the ccall estimate as estimated by SD/B

   - ccgroup: number of correctly classified samples from each group

   - ngroup: number of samples from each group

*/





  int i, j, *dpred, *d, *groupspred, ncolpred, ncolsumx, zero=0, one=1, *nobsgroup, *cluslist;

  double *xpred, *apred, *lpred, *vpred, *sumxpred, *prodxpred, *nobsxpred, cv, lhood, pcum, ccobs, seccobs, *ccgroupobs;



  ncolpred= (*deltat)*(*K);

  xpred= dvector(0,(*nrow)*ncolpred);

  dpred= ivector(0,(*nrow));

  d= ivector(0,(*nrow));

  if (*equalcv == 1) apred= dvector(0,*nrow); else apred= dvector(0,(*nrow)*(*K));

  lpred= dvector(0,(*nrow)*(*K));

  vpred= dvector(0,(*nrow)*(*npat));

  ccgroupobs= dvector(0,*K); nobsgroup= ivector(0,*K);

  cluslist= ivector(0,*nclust);



  for (i=0; i<(*nclust); i++) { cluslist[i]= i; }

  cluslist[*nclust]= -1;



  groupspred= ivector(0,ncolpred);                                            //indicate group each new obs will belong to

  for (i=0; i<ncolpred; i++) { groupspred[i]= floor(i / (*deltat)); }



  if (*usesel == 0) { sel_mostDEgenes(nsel,sel,genelimit,v0thre,nrow,npat,v); }  //select most DE genes



  ncolsumx= 0; for (i=0; i<(*npat); i++) { ncolsumx += ngrouppat[i]; }

  if ((*usesumx)==0) { compute_sumxC(sumx,prodx,nobsx,equalcv,nsel,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one); }

  sumxpred= dvector(0,(*nrow)*ncolsumx); nobsxpred= dvector(0,ncolsumx);

  if (*equalcv) {

    prodxpred= dvector(0,(*nrow));

  } else {

    prodxpred= dvector(0,(*nrow)*ncolsumx);;

  }



  (*ccall)= (*seccall)= 0; for (i=0; i<(*K); i++) { ccgroup[i]= 0; ngroup[i]= 0; }

  cv= 1; pcum= Kprob[0];



  for (i=0; i<(*B) && (cv>(*preceps) || i<100); i++) {                                        //loop over simulations

    simpred_ggC(xpred,dpred,apred,lpred,&zero,deltat,groups,K,nsel,sel,nrow,ncol,x,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,v,npat,patterns,ngrouppat,sumx,prodx,nobsx,&one,gapprox);

    compute_sumxC(sumxpred,prodxpred,nobsxpred,equalcv,nsel,sel,&ncolsumx,&ncolpred,xpred,groupspred,K,npat,patterns,ngrouppat,&one); //suff stat for xpred

    pp_ggC(vpred,&lhood,nsel,sel,&ncolpred,xpred,groupspred,K,alpha0,nu,balpha,nualpha,equalcv,nclust,cluslist,rho,prob,npat,patterns,ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,&one,gapprox); //post prob of patterns

    add_sumxC(sumxpred,prodxpred,nobsxpred,equalcv,nsel,sel,&ncolsumx,sumx,prodx,nobsx); //add sumx and sumxpred into sumxpred



    utsample_ggC(&ccobs,&seccobs,ccgroupobs,nobsgroup,&one,preceps,genelimit,v0thre,nsel,sel,&one,nrow,&ncolpred,xpred,groupspred,vpred,K,Kprob,alpha0,nu,balpha,nualpha,equalcv,nclust,rho,prob,npat,patterns,ngrouppat,&ncolsumx,sumxpred,prodxpred,nobsxpred,&one,gapprox); //observed utility (only 1 sample taken)



    (*ccall) += ccobs; for (j=0; j<(*K); j++) { ccgroup[j] += ccgroupobs[j]; ngroup[j] += nobsgroup[j]; }

    (*seccall) += ccobs*ccobs;

    if (i>=50) { cv= sqrt((*seccall)/(i+1) - (*ccall)*(*ccall)/((i+1)*(i+1)))/sqrt(i+i); }

  }                                                                                            //end loop over simulations

  (*ccall)= (*ccall) / i;

  (*seccall)= sqrt((*seccall)/i - (*ccall)*(*ccall))/sqrt(0.0+i); 



  free_ivector(cluslist,0,*nclust);

  free_dvector(xpred,0,(*nrow)*ncolpred);

  free_ivector(dpred,0,(*nrow));

  free_ivector(d,0,(*nrow));

  if (*equalcv == 1) free_dvector(apred,0,*nrow); else free_dvector(apred,0,(*nrow)*(*K));

  free_dvector(lpred,0,(*nrow)*(*K));

  free_dvector(vpred,0,(*nrow)*(*npat));

  free_ivector(groupspred,0,ncolpred);

  free_dvector(sumxpred,0,(*nrow)*ncolsumx);

  if (*equalcv == 1) {

    free_dvector(prodxpred,0,(*nrow));

  } else {

    free_dvector(prodxpred,0,(*nrow)*ncolsumx);

  }

  free_dvector(nobsxpred,0,ncolsumx);

  free_dvector(ccgroupobs,0,*K); free_ivector(nobsgroup,0,*K);



}







/*********************************************************************************************

                                    POSTERIOR COMPUTATIONS

*********************************************************************************************/



void compute_sumxC(double *sumx, double *prodx, double *nobsx, int *equalcv, int *nsel, int *sel, int *ncolsumx, int *ncol, double *x, int *groups, int *K, int *npat, int *patterns, int *ngrouppat, int *init0) {

/* Computes sufficient statistics for Gamma/Gamma model: sums of columns of x and number of terms in the sum */

/* Input arguments

  - equalcv: set to 1 to indicate constant CV across groups, 0 otherwise

  - nsel: sufficient statistics will only be computed for nsel genes

  - sel: vector with its first nsel positions indicating the genes to compute sufficient statistics for

  - ncolsumx: number of columns of sumx

  - ncol: number of columns of x

  - x: matrix with observations (genes in rows, samples in columns).

  - groups: vector indicating what group each column in x corresponds to

  - K: number of groups e.g. group 0 for control, group 1 for type 1 cancer, group 2 for type B cancer

  - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

  - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

  - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

               (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as argument)

  - init0: if *init0==1 sumx, prodx and nobsx are initialized to 0, otherwise the corresponding sums are added to its previous contents

*/

/* Output arguments

  - sumx: sums of columns of x for each gene in sel and each group in each expression pattern (nrow rows, ncolsumx cols). The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on. Only rows sel[0],sel[1] etc are computed, the rest are left untouched.

  - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

  - nobsx: vector with number of terms in the sum. The first ngrouppat[0] elements correspond to pattern 0 and so on.

*/



  int i, j, k, di, colini, group;



  //Compute sums

  if (*equalcv == 1) {

    for (i=0; i<(*nsel); i++) {                                         //for each gene

      if ((*init0)==1) {

	for (k=0; k<(*ncolsumx); k++) { sumx[sel[i]*(*ncolsumx)+k] = 0; }

	prodx[sel[i]] = 0;

      }

      for (j=0; j<(*ncol); j++) {                                       //for each sample

	for (di=0,colini=0; di<(*npat); di++) {                                        //for each pattern

	  group= patterns[di*(*K)+groups[j]];                                //find cluster of groups each obs corresponds to

	  if (x[sel[i]*(*ncol)+j]>0.0001) { sumx[sel[i]*(*ncolsumx)+colini+group] += x[sel[i]*(*ncol)+j]; } else { sumx[sel[i]*(*ncolsumx)+colini+group] += 0.0001; }

	  colini += ngrouppat[di];

	}

	if (x[sel[i]*(*ncol)+j]>0.0001) { prodx[sel[i]] += log(x[sel[i]*(*ncol)+j]); } else { prodx[sel[i]] += log(0.0001); }

      }  

    }

  } else {

    for (di=0,colini=0; di<(*npat); di++) {                                        //for each pattern

      for (i=0; i<(*nsel); i++) {                                         //for each gene

	if ((*init0)==1) { 

	  for (j=0; j<ngrouppat[di]; j++) { sumx[sel[i]*(*ncolsumx)+colini+j] = prodx[sel[i]*(*ncolsumx)+colini+j] = 0; }

	}

	for (j=0; j<(*ncol); j++) {                                       //for each sample

	  group= patterns[di*(*K)+groups[j]];                                //find cluster of groups each obs corresponds to

	  if (x[sel[i]*(*ncol)+j]>0.0001) {

	    sumx[sel[i]*(*ncolsumx)+colini+group] += x[sel[i]*(*ncol)+j];

	    prodx[sel[i]*(*ncolsumx)+colini+group] += log(x[sel[i]*(*ncol)+j]);

	  } else {

	    sumx[sel[i]*(*ncolsumx)+colini+group] += 0.0001;

	    prodx[sel[i]*(*ncolsumx)+colini+group] += log(0.0001);

	  }

	}                                                                 //end for each sample

      }                                                                   //end for each gene

      colini += ngrouppat[di];

    }

  }                                                                     //end for each pattern



  //Compute number of terms in the sums

  colini= 0;

  for (di=0; di<(*npat); di++) {                                        //for each pattern

      if ((*init0)==1) { 

	for (j=0; j<ngrouppat[di]; j++) { nobsx[colini+j] = 0; }

      }

      for (j=0; j<(*ncol); j++) {                                       //for each sample

	group= patterns[di*(*K)+groups[j]];                                //find cluster of groups each obs corresponds to

	nobsx[colini+group] += 1;

      }                                                                 //end for each sample

    colini += ngrouppat[di];

  }                                                                     //end for each pattern



}



void add_sumxC(double *sumxnew, double *prodxnew, double *nobsnew, int *equalcv, int *nsel, int *sel, int *ncolsumx, double *sumx, double *prodx, double *nobsx) {

/* Adds contents of sumx, prodx, nobsx to sumxnew, prodxnew, nobsnew */

/* Input arguments

   - equalcv: set to 1 to indicate constant CV across groups, 0 otherwise

   - nsel: contents will only be copied for nsel genes

   - sel: vector with its first nsel positions indicating the genes to make copies for

   - ncol: number of columns of sumx, and also number of elements in nobsx

   - sumx: matrix to be copied from

   - prodx: matrix to be copied from

   - nobsx: vector to be copied from

*/

/* Ouput arguments

   - sumxnew, prodxnew, nobsnew: copies of sumx, prodx and nobsx. Only the nsel first positions indicated in sel are actually copied.

*/



  int i,j;



  if (*equalcv == 1) {

    for (i=0; i<(*nsel); i++) {

      for (j=0; j<(*ncolsumx); j++) { sumxnew[sel[i]*(*ncolsumx)+j]+= sumx[sel[i]*(*ncolsumx)+j]; }

      prodxnew[sel[i]]+= prodx[sel[i]];

    }

  } else {

    for (i=0; i<(*nsel); i++) {

      for (j=0; j<(*ncolsumx); j++) {

	sumxnew[sel[i]*(*ncolsumx)+j]+= sumx[sel[i]*(*ncolsumx)+j]; prodxnew[sel[i]*(*ncolsumx)+j]+= prodx[sel[i]*(*ncolsumx)+j];

      }

    }

  }

  for (j=0; j<(*ncolsumx); j++) { nobsnew[j]+= nobsx[j]; }

}



void copy_sumxC(double *sumxnew, double *prodxnew, double *nobsnew, int *equalcv, int *nsel, int *sel, int *ncolsumx, double *sumx, double *prodx, double *nobsx) {

/* Copies contents of sumx, prodx, nobsx into sumxnew, prodxnew, nobsnew */

/* Input arguments

   - equalcv: set to 1 to indicate constant CV across groups, 0 otherwise

   - nsel: contents will only be copied for nsel genes

   - sel: vector with its first nsel positions indicating the genes to make copies for

   - ncol: number of columns of sumx, and also number of elements in nobsx

   - sumx: matrix to be copied from

   - prodx: matrix to be copied from

   - nobsx: vector to be copied from

*/

/* Ouput arguments

   - sumxnew, prodxnew, nobsnew: copies of sumx, prodx and nobsx. Only the nsel first positions indicated in sel are actually copied.

*/



  int i,j;



  if (*equalcv == 1) {

    for (i=0; i<(*nsel); i++) { for (j=0; j<(*ncolsumx); j++) { sumxnew[sel[i]*(*ncolsumx)+j]= sumx[sel[i]*(*ncolsumx)+j]; } }

    prodxnew[sel[i]]= prodx[sel[i]];

  } else {

    for (i=0; i<(*nsel); i++) { 

      for (j=0; j<(*ncolsumx); j++) {

	sumxnew[sel[i]*(*ncolsumx)+j]= sumx[sel[i]*(*ncolsumx)+j]; prodxnew[sel[i]*(*ncolsumx)+j]= prodx[sel[i]*(*ncolsumx)+j];

      }

    }

  }

  for (j=0; j<(*ncolsumx); j++) { nobsnew[j]= nobsx[j]; }



}





void pp_ggC(double *v, double *lhood, int *nsel, int *sel, int *ncol, double *x, int *groups, int *K, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, int *cluslist, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, double *sumxpred, double *prodxpred, double *nobsxpred, int *usesumx, int *gapprox) {

/* Computes posterior probabilities of each expression pattern given in 'patterns' conditional on data x and hyperparameters for Gamma/Gamma model */

/*   Input arguments

   - nsel: posterior probabilities are computed for nsel genes only.

   - sel: vector containing the keys of the nsel genes

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: matrix with observations used to fit the model (doesn't need to be in any specific order)

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - alpha0,nu: prior for lambda is Gamma(alpha0,nu)

   - balpha,nualpha: prior for alpha is Gamma(balpha,nualpha)

   - equalcv: set to 1 to indicate constant CV across groups, 0 otherwise

   - nclust: number of clusters in hyper-prior for (alpha,lambda)

   - cluslist: list of clusters to be taken into account when computing post prob. This option can be used to exclude clusters with very small probabilities

   - rho: mixing probabilities for clusters in hyper-prior for (alpha,lambda)

   - prob: vector with estimated probability of each expression pattern (length *npat)

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on. If usesumx==0 its values are computed but it still must be a vector of the appropriate size.

   - prodx: same as sumx but contains logarithm of the products i.e. sums of logarithms

   - nobsx: number of terms in the sum. As for sumx, the first ngrouppat[0] cols correspond to pattern 0 and so on.

   - sumxpred: same as sumx but can be used to store sums for obs from the predictive.

   - prodxpred: same as prodx but can be used to store log prods for obs from the predictive.

   - nobsxpred: same as nobsx but can be used to store number of terms in sums from the predictive.

   - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, prodx, nobsx, sumxpred, prodxpred, nobsxpred and ignore the argument x (actually only the sums sumx+sumxpred, prodx+prodxpred and nobs+nobsxpred are used). For usesumx==0 they're computed from the data x and sumx, prodx, nobsx, sumxpred, prodxpred, nobsxpred are ignored.

   - gapprox: for gapprox==1 the normalization constant for the conjugate Gamma shape distributions is computed approximately



   Output arguments

   - v: matrix with post prob of each expression pattern (*nrow rows,*npat cols). In 2 group case reduces to post prob of equal and differential expression

   - lhood: log-likelihood of the observed data given the hyper-parameter values

*/



  int i, j, k, m, *colini, ncolsumx, init0=1, usexpred=1, firstclus;

  double vsum, r0=0, rcur, rsum, lgene;

  //  FILE *ifile; //debug

  //ifile= openOut("lhood.txt"); //debug



colini= ivector(0,*npat);

colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];



if ((*usesumx==0)) {                                               //if suff stat not pre-computed

  compute_sumxC(sumx,prodx,nobsx,equalcv,nsel,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&init0);

  for (j=0; j<ncolsumx; j++) {

    for (i=0; i<(*nsel); i++) { sumxpred[i*ncolsumx+j]= 0; }

    nobsxpred[j]= 0;

  }

  if (*equalcv==1) {

    for (i=0; i<(*nsel); i++) { prodxpred[i]= 0; }

  } else {

    for (j=0; j<ncolsumx; j++) { for (i=0; i<(*nsel); i++) { prodxpred[i*ncolsumx+j]= 0; } }

  }

}



(*lhood)= 0;

for (i=0; i<(*nsel); i++) {                                      //for each gene

  lgene= 0;

  for (j=0, vsum=1; j<(*npat); j++) {                            //for each expression pattern

    for (k=0, rsum=1, firstclus=1; cluslist[k]!=(-1); k++) {     //for each cluster

      m= cluslist[k];

      rcur= pdfcond_pat_clus(sel[i],j,m,alpha0,nu,balpha,nualpha,ngrouppat,colini,ncolsumx,sumx,sumxpred,prodx,prodxpred,nobsx,nobsxpred,equalcv,usexpred);

      //      if (j==(*npat -1)) fprintf(ifile,"%.6f \n",rcur); else fprintf(ifile,"%.6f ",rcur); //debug

      lgene += exp(rcur)*prob[j]*rho[m];

      if (firstclus==1) { r0= rcur; firstclus= 0; } else { rsum += exp(rcur-r0)*rho[m]/rho[0]; }

    }

    v[sel[i]*(*npat)+j]= log(rsum) + r0 + log(rho[0]);

    if (j>0) { 

      v[sel[i]*(*npat)+j]= min_xy(exp(v[sel[i]*(*npat)+j] - v[sel[i]*(*npat)])*prob[j]/prob[0],exp(700.0)); //cut at exp(700) to avoid overflow

      vsum+= v[sel[i]*(*npat)+j];

    }

  }

  v[sel[i]*(*npat)]= 1/vsum;                                     //normalize post probabilities

  for (j=1; j<(*npat); j++) { v[sel[i]*(*npat)+j]= v[sel[i]*(*npat)+j] / vsum; }

  (*lhood) += log(lgene);

  //  fprintf(ifile,"%.6f \n",log(lgene)); //debug

}                                                                //end for each gene



//fclose(ifile); //debug

free_ivector(colini,0,*npat);



}





double pdfcond_pat_clus(int geneid, int patid, int clusid, double *alpha0, double *nu, double *balpha, double *nualpha, int *ngrouppat, int *colini, int ncolsumx, double *sumx, double *sumxpred, double *prodx, double *prodxpred, double *nobsx, double *nobsxpred, int *equalcv, int usexpred) {

//Density of data for gene geneid conditional on it following expression pattern patid and belonging to hyper-parameter cluster clusid

  int l, newton= 2, logscale=1, lengtha;

  double rcur, *n, *s, pr, b1, b2;



  if (*equalcv == 1) {

    rcur= ngrouppat[patid]*(alpha0[clusid]*log(alpha0[clusid]/nu[clusid]) -gamln(alpha0+clusid)) + balpha[0]*log(balpha[0]/nualpha[0]) -gamln(balpha);

    lengtha= ngrouppat[patid];

    pr= prodx[geneid];

    if (usexpred==0) {

      n= nobsx+colini[patid];

      s= sumx+geneid*ncolsumx+colini[patid];

      b1= balpha[0]/nualpha[0] - pr; b2= alpha0[clusid]/nu[clusid];

      rcur += kcgammaC(n,balpha,&b1,alpha0+clusid,&b2,s,&lengtha,&newton,&logscale);

    } else {

      n= dvector(0,ngrouppat[patid]); s= dvector(0,ngrouppat[patid]);

      for (l=0; l<ngrouppat[patid]; l++) { 

	n[l]= nobsx[colini[patid]+l] + nobsxpred[colini[patid]+l];

	s[l]= sumx[geneid*ncolsumx+colini[patid]+l] + sumxpred[geneid*ncolsumx+colini[patid]+l];

      }

      b1= balpha[0]/nualpha[0] - pr; b2= alpha0[clusid]/nu[clusid];

      rcur += kcgammaC(n,balpha,&b1,alpha0+clusid,&b2,s,&lengtha,&newton,&logscale);

      free_dvector(n,0,ngrouppat[patid]); free_dvector(s,0,ngrouppat[patid]);

    }

  } else {

    n= dvector(0,1); s= dvector(0,1);

    lengtha= 1;

    rcur= ngrouppat[patid]*(alpha0[clusid]*log(alpha0[clusid]/nu[clusid]) -gamln(alpha0+clusid) + balpha[0]*log(balpha[0]/nualpha[0]) -gamln(balpha));

    for (l=0; l<(ngrouppat[patid]); l++) {

      n[0]= nobsx[colini[patid]+l];

      s[0]= sumx[geneid*ncolsumx+colini[patid]+l];

      pr= prodx[geneid*ncolsumx+colini[patid]+l];

      if (usexpred==1) {

	n[0] += nobsxpred[colini[patid]+l];

	s[0] += sumxpred[geneid*ncolsumx+colini[patid]+l];

	pr += prodxpred[geneid*ncolsumx+colini[patid]+l];

      }

      b1= balpha[0]/nualpha[0] - pr; b2= alpha0[clusid]/nu[clusid];

      //if (geneid==0) printf("Gene %d, Pattern %d. %.4f %.4f %.4f %.4f %.4f %.4f \n",geneid,patid,n[0],balpha[0],b1,alpha0[clusid],b2,s[0]); //debug

      rcur += kcgammaC(n,balpha,&b1,alpha0+clusid,&b2,s,&lengtha,&newton,&logscale);

    }

    free_dvector(n,0,1); free_dvector(s,0,1);

  }



  return(rcur);

}





void estep(double *probest, double *rhoest, double *lhood, int *justlhood, int *nsel, int *sel, int *K, double *logalpha0, double *lognu, double *logbalpha, double *lognualpha, int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, int *colini, int *ncolsumx, double *sumx, double *prodx, double *nobsx, int *gapprox) {

//E-step within the EM algorithm. Computes expected value of mixing probabilities and cluster probabilities, given the hyper-parameters and the data

/* Input: all as in pp_ggC, except for the following:

   - justlhood: if justlhood==1, marginal likelihood is evaluated but probest and rhoest are not calculated.

   - logalpha0,lognu,logbalpha,lognualpha: as in pp_ggC, except that the parameters are in log-scale

   - logitrho, logitprob: as in pp_ggC, except that it's in logit scale



   Output

   - probest[0..*npat-1]: expected value of expression pattern probabilities

   - rhoest[0..*nclust-1]: expected value of cluster probabilities

   - lhood: marginal likelihood

*/



  int i, j, k;

  double *lgene, ***r, **patclus, sumprobest, sumrhoest, *alpha0, *nu, *balpha, *nualpha;

  //  FILE *ifile; //debug



  //put parameter values back to original scale

  alpha0= dvector(0,*nclust); nu= dvector(0,*nclust);

  balpha= dvector(0,1); nualpha= dvector(0,1);



  for (i=0; i< *nclust; i++) { alpha0[i]= exp(logalpha0[i]); nu[i]= exp(lognu[i]); }

  balpha[0]= exp(logbalpha[0]); nualpha[0]= exp(lognualpha[0]);   



  //allocate memory

  r= darray3(*nsel,*npat,*nclust);

  lgene= dvector(0,*nsel);

  patclus= dmatrix(0,*npat-1,0,*nclust-1);



  //Evaluate log-likelihood

  (*lhood)= 0;

  for (i=0; i<(*nsel); i++) {                                      //for each gene

    for (j=0, lgene[i]=0; j<(*npat); j++) {                        //for each expression pattern

      for (k=0; k< *nclust; k++) {                                 //for each cluster

	//        if (__iglobal==38) printf("i=%d, j=%d \n",i,j); //debug

	r[sel[i]][j][k]= exp(pdfcond_pat_clus_nopred(sel[i],j,k,alpha0,nu,balpha,nualpha,ngrouppat,colini,*ncolsumx,sumx,prodx,nobsx,equalcv))*prob[j]*rho[k];

	lgene[i] += r[sel[i]][j][k] + 1.0E-30;  //avoid 0 likelihood by adding offset

      }

    }

    (*lhood) += log(lgene[i]);

  }



  //ifile= openOut("lhood.dat"); //debug

  //for (i=0; i<(*nsel); i++) { fprintf(ifile,"%d %.8f \n",i,lgene[i]); } //debug

  //fclose(ifile); //debug



  //E-step

  if (*justlhood ==0) {

    //sum posterior probabilities of (pattern,cluster) across genes

    sumprobest= sumrhoest= 0;

    for (j=0; j<(*npat); j++) {         

      for (k=0; k< *nclust; k++) {      

        for (i=0, patclus[j][k]=0; i<(*nsel); i++) {

	  patclus[j][k]+= r[sel[i]][j][k]/lgene[sel[i]];

	}

	probest[j]+= patclus[j][k]; sumprobest+= patclus[j][k];

	rhoest[k]+= patclus[j][k]; sumrhoest+= patclus[j][k];

      }

    }

    //normalize

    for (j=0; j<(*npat); j++) probest[j]/= sumprobest;

    for (k=0; k< *nclust; k++) rhoest[k]/= sumrhoest;

  }



  free_dvector(alpha0,0,*nclust); free_dvector(nu,0,*nclust);

  free_dvector(balpha,0,1); free_dvector(nualpha,0,1);



  free_darray3(r,*nsel,*npat,*nclust);

  free_dvector(lgene,0,*nsel);

  free_dmatrix(patclus,0,*npat-1,0,*nclust-1);



}







void lhoodnopat(double *probest, double *rhoest, double *lhood, int *justlhood, int *nsel, int *sel, int *K, double *logalpha0, double *lognu, double *logbalpha, double *lognualpha, int *equalcv, int *nclust, double *rho, double *prob, int *npat, int *patterns, int *ngrouppat, int *colini, int *ncolsumx, double *sumx, double *prodx, double *nobsx, int *gapprox) {

//log-likelihood under the complete null hypothesis that all genes are EE

/* Input: all as in estep



   Output

   - lhood: marginal likelihood

*/



  int i, j, k;

  double *lgene, ***r, **patclus, *alpha0, *nu, *balpha, *nualpha;

  //  FILE *ifile; //debug



  //put parameter values back to original scale

  alpha0= dvector(0,*nclust); nu= dvector(0,*nclust);

  balpha= dvector(0,1); nualpha= dvector(0,1);



  for (i=0; i< *nclust; i++) { alpha0[i]= exp(logalpha0[i]); nu[i]= exp(lognu[i]); }

  balpha[0]= exp(logbalpha[0]); nualpha[0]= exp(lognualpha[0]);   



  //allocate memory

  r= darray3(*nsel,*npat,*nclust);

  lgene= dvector(0,*nsel);

  patclus= dmatrix(0,*npat-1,0,*nclust-1);



  //Evaluate log-likelihood

  (*lhood)= 0;

  for (i=0; i<(*nsel); i++) {                                      //for each gene

    for (j=0, lgene[i]=0; j<1; j++) {                        //only for null expression pattern

      for (k=0; k< *nclust; k++) {                                 //for each cluster

	r[sel[i]][j][k]= exp(pdfcond_pat_clus_nopred(sel[i],0,k,alpha0,nu,balpha,nualpha,ngrouppat,colini,*ncolsumx,sumx,prodx,nobsx,equalcv))*prob[j]*rho[k];

	lgene[i] += r[sel[i]][j][k] + 1.0E-30;  //avoid 0 likelihood by adding offset

      }

    }

    (*lhood) += log(lgene[i]);

  }



  free_dvector(alpha0,0,*nclust); free_dvector(nu,0,*nclust);

  free_dvector(balpha,0,1); free_dvector(nualpha,0,1);



  free_darray3(r,*nsel,*npat,*nclust);

  free_dvector(lgene,0,*nsel);

  free_dmatrix(patclus,0,*npat-1,0,*nclust-1);



}









double pdfcond_pat_clus_nopred(int geneid, int patid, int clusid, double *alpha0, double *nu, double *balpha, double *nualpha, int *ngrouppat, int *colini, int ncolsumx, double *sumx, double *prodx, double *nobsx, int *equalcv) {

//Same as pdfcond_pat, but runs faster as it doesn't consider data from the predictive sumxpred, prodxpred and nobsxpred.

  int l, newton= 2, logscale=1, lengtha;

  double rcur, *n, *s, pr, b1, b2;



  if (*equalcv == 1) {

    //    printf("2: %.5f \n",balpha[0]*log(balpha[0]/nualpha[0]));

    //printf("3: %.5f \n",-gamln(balpha));

    //printf("1: %.5f \n",ngrouppat[patid]*(alpha0[clusid]*log(alpha0[clusid]/nu[clusid]) -gamln(alpha0+clusid)));



    rcur= ngrouppat[patid]*(alpha0[clusid]*log(alpha0[clusid]/nu[clusid]) -gamln(alpha0+clusid)) + balpha[0]*log(balpha[0]/nualpha[0]) -gamln(balpha);

    lengtha= ngrouppat[patid];

    pr= prodx[geneid];

    n= nobsx+colini[patid];

    s= sumx+geneid*ncolsumx+colini[patid];

    b1= balpha[0]/nualpha[0] - pr; b2= alpha0[clusid]/nu[clusid];

    rcur += kcgammaC(n,balpha,&b1,alpha0+clusid,&b2,s,&lengtha,&newton,&logscale);

  } else {

    n= dvector(0,1); s= dvector(0,1);

    lengtha= 1;

    rcur= ngrouppat[patid]*(alpha0[clusid]*log(alpha0[clusid]/nu[clusid]) -gamln(alpha0+clusid) + balpha[0]*log(balpha[0]/nualpha[0]) -gamln(balpha));

    for (l=0; l<(ngrouppat[patid]); l++) {

      n[0]= nobsx[colini[patid]+l];

      s[0]= sumx[geneid*ncolsumx+colini[patid]+l];

      pr= prodx[geneid*ncolsumx+colini[patid]+l];

      b1= balpha[0]/nualpha[0] - pr; b2= alpha0[clusid]/nu[clusid];

      rcur += kcgammaC(n,balpha,&b1,alpha0+clusid,&b2,s,&lengtha,&newton,&logscale);

    }

    free_dvector(n,0,1); free_dvector(s,0,1);

  }



  return(rcur);

}





void posmeans_ggC(double *posmeans, int *underpattern, int *K, int *nsel, int *sel, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *gapprox) {

/* Posterior group means. Approximation based on setting the shape parameter to its posterior mean.

   Input arguments

   - underpattern: expression patterns assumed to be true. Posterior mean are computed under this pattern.

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - nsel: expression values are simulated for nsel genes only.

   - sel: vector containing the keys of the nsel genes

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - balpha: estimate for b parameter in hyper-prior for alpha

   - v: matrix with post prob of each expression pattern (*nrow rows,*npat cols). In 2 group case reduces to post prob of equal and differential expression

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

  - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on.

  - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

  - nobsx: vector with number of terms in the sum. The first ngrouppat[0] elements correspond to pattern 0 and so on.

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used

   Output arguments

   - posmeans: matrix of nsel rows and K cols with posterior group means.

*/



  int i, k, m, di, ncolsumx, *colini, coldi, newton=2, usexpred=0, lengtha;

  double rsum, a, b1, b2, s, *vclus, *sumxpred, *prodxpred, *nobsxpred, aest, best;



  sumxpred= dvector(0,1); prodxpred= dvector(0,1); nobsxpred= dvector(0,1); //never used: initialized to avoid warnings



  vclus= dvector(0,*nclust);

  colini= ivector(0,*npat);

  colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

  ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];



  di= (*underpattern);

  for (i=0; i<(*nsel); i++) {                                 //compute post means for nsel genes only

    for (k=0; k<(*K); k++) {

      //cluster probabilities

      if (*nclust>1) {

        for (m=0, rsum=1; m<(*nclust); m++) {

  	  vclus[m]= pdfcond_pat_clus(i,di,m,alpha0,nu,balpha,nualpha,ngrouppat,colini,ncolsumx,sumx,sumxpred,prodx,prodxpred,nobsx,nobsxpred,equalcv,usexpred) + log(rho[m]);

	  if (m>0) { vclus[m]= exp(vclus[m]-vclus[0]); rsum+= vclus[m]; }

        }

        vclus[0]= 1/rsum;

        for (m=1; m<(*nclust); m++) { vclus[m]= vclus[m]/rsum; }

      } else {

        vclus[0]= 1.0;

      }



      //estimate alpha,l 

      if (*equalcv == 1) {

        lengtha= ngrouppat[di];

        coldi= colini[di];

        b1= balpha[0]/nualpha[0] - prodx[sel[i]];

        for (m=0, posmeans[i*(*K)+k]=0; m<(*nclust); m++) {

	  b2= alpha0[m]/nu[m];

	  gapprox_par(&aest,&best,nobsx+coldi,balpha,&b1,alpha0+m,&b2,sumx+sel[i]*ncolsumx+coldi,&lengtha,&newton);

	  a= aest/best;

	  aest= alpha0[m] + a*nobsx[coldi+patterns[di*(*K)+k]];

	  best= alpha0[m]/nu[m] + a*sumx[sel[i]*ncolsumx+coldi+patterns[di*(*K)+k]];

	  if (aest>1) posmeans[i*(*K)+k]+= vclus[m]*best/(aest-1.0); else posmeans[i*(*K)+k]+= vclus[m]*best/(aest+1.0);

        }

      } else {

        lengtha= 1;

        coldi= colini[di]+patterns[di*(*K)+k];

        s= sumx[sel[i]*ncolsumx + coldi];

        b1= balpha[0]/nualpha[0] - prodx[sel[i]*ncolsumx + coldi];

        for (m=0; m<(*nclust); m++) { 

	  b2= alpha0[m]/nu[m];

	  gapprox_par(&aest,&best,nobsx+coldi,balpha,&b1,alpha0+m,&b2,&s,&lengtha,&newton);

          a= aest/best;

	  aest= alpha0[m]/nu[m] + a*s;

	  best= alpha0[m] + a*nobsx[coldi];

	  if (aest>1) posmeans[i*(*K)+k]+= vclus[m]*best/(aest-1.0); else posmeans[i*(*K)+k]+= vclus[m]*best/(aest+1.0);

        }

      }

    } //End k for

  }  //End i for



  free_ivector(colini,0,*npat);

  free_dvector(vclus,0,*nclust);

  free_dvector(sumxpred,0,1); free_dvector(prodxpred,0,1); free_dvector(nobsxpred,0,1);



}





/*********************************************************************************************

                                  SIMULATE FROM THE PREDICTIVE

*********************************************************************************************/



void simhyperpar_ggC(double *alpha0, double *nu, double *balpha, double *nualpha, int *nclust, double *rho, double *prob, int *cluslist, int *equalcv, double *a_alpha0, double *b_alpha0, double *a_nu, double *b_nu, double *a_balpha, double *b_balpha, double *a_nualpha, double *b_nualpha, double *p_rho, double *p_prob, int *nrow, double *sumd, double *sumci, double *ngroupstot, double *sumalpha, double *sumlogalpha, double *suminvlambda, double *sumlambda, double *sumloglambda, int *npat, int *ngrouppat, int *gapprox, int *settomean) {

/* Simulates hyper-parameter values from the posterior conditional on all other parameters in the Gamma/Gamma model */

/* Input arguments

   - equalcv:

   - a_alpha0, b_alpha0: prior for alpha0 is Gamma(a_alpha0,b_alpha0)

   - a_nu, b_nu: prior for nu is Gamma(a_nu,b_nu)

   - a_balpha, b_balpha: prior for balpha is Gamma(a_balpha,b_balpha)

   - a_nualpha, b_nualpha: prior for nualpha is Gamma(a_nualpha,b_nualpha)

   - p_rho: prior for rho is Dirichlet(p_rho). p_rho is a vector of length nclust

   - p_prob: prior for prob is Dirichlet(p_prob). p_prob is a vector of length npat

   - nrow: number of genes (length of vector d)

   - sumd: vector of length *npat indicating the number of genes falling into each expression pattern.

   - ngroupstot: number of distinct (alpha,lambda) parameters

   - sumlambda: sum of distinct lambda parameters for all genes

   - sumloglambda: sum of distinct log(lambda) parameters

   - sumalpha: sum of distinct alpha parameters for all genes

   - sumlogalpha: sum of distinct log(alpha) parameters for all genes

   - suminvlambda: sum of 1/lambda

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used

   - settomean: for settomean==1 hyper-par values are set equal to their expected value (condit on all other params)

   Output arguments

   - alpha0: value drawn for alpha0

   - nu: value drawn for nu

   - balpha: value drawn for balpha

   - nualpha: value drawn for nualpha

   - prob: value drawn for prob (mixing probabilities)

   - cluslist: list of clusters that have at least one observation assigned to them

*/



  int i, m, one=1, newton=5, lengtha=1;

  double b1, sumla, ntot, suma, sumtot, aest, best;



  sumla= ntot= suma= 0;

  for (i=0, m=0; i<(*nclust); i++) {

    if (sumci[i]>0) {

      cluslist[m]= i; m++;

      b1= *b_alpha0 + sumloglambda[i];

      rcgammaC(alpha0+i,&one,ngroupstot+i,a_alpha0,&b1,a_nu,b_nu,suminvlambda+i,&lengtha,&newton); //draw alpha0

      nu[i]= 1.0/gengam(*b_nu + alpha0[i]*suminvlambda[i],*a_nu + alpha0[i]*ngroupstot[i]); //draw nu

    } else {

      //alpha0[i]= (*a_alpha0)/(*b_alpha0);  //expectation a priori

      //nu[i]= (*b_nu)/(*a_nu + 1.0); //mode a priori

      alpha0[i]= gengam(*b_alpha0,*a_alpha0);  //generate from prior

      nu[i]= 1.0/gengam(*b_nu,*a_nu);

    }

    if (alpha0[i]<=.0001) alpha0[i]= (*a_alpha0)/(*b_alpha0);

    if (nu[i]>exp(500)) nu[i]= (*b_nu)/(*a_nu + 1.0);

    sumla += sumlogalpha[i];

    ntot += ngroupstot[i];

    suma += sumalpha[i];

  }

  cluslist[m]= -1;

  b1= *b_balpha - sumla;

  if (*equalcv == 1) { ntot= *nrow; }



  if (*settomean==1) {

    gapprox_par(&aest,&best,&ntot,a_balpha,&b1,a_nualpha,b_nualpha,&suma,&lengtha,&newton); 

    (*balpha)= aest/best;  //expected value of balpha

    if (*a_nualpha+balpha[0]*ntot>1) (*nualpha)= (*b_nualpha+balpha[0]*suma)/(*a_nualpha+balpha[0]*ntot-1.0); else (*nualpha)= (*b_nualpha+balpha[0]*suma)/(*a_nualpha+balpha[0]*ntot+1.0); //expected value of nualpha. If infinite, set to mode.

    if (*nclust > 1) {

      for (i=0, sumtot=0; i<*nclust; i++) { sumci[i]+= p_rho[i]; sumtot+=sumci[i]; } 

      for (i=0; i<*nclust; i++) rho[i]= sumci[i]/sumtot; //expected value of cluster prob

    } else {

      rho[0]= 1;

    }

    for (i=0, sumtot=0; i<*npat; i++) { sumd[i]+= p_prob[i]; sumtot+= sumd[i]; } 

    for (i=0; i<*npat; i++) { prob[i]= sumd[i]/sumtot; } //expected value of mixing prob

  } else {   //settomean==0

    rcgammaC(balpha,&one,&ntot,a_balpha,&b1,a_nualpha,b_nualpha,&suma,&lengtha,&newton); //draw balpha

    nualpha[0]= 1.0/gengam(*b_nualpha+balpha[0]*suma,*a_nualpha+balpha[0]*ntot); //draw nualpha

    if (*nclust > 1) {

      for (i=0; i<*nclust; i++) { sumci[i]+= p_rho[i]; }                      //draw cluster prob values

      rdirichlet(rho,sumci,nclust);

    } else {

      rho[0]= 1;

    }

    for (i=0; i<*npat; i++) { sumd[i]+= p_prob[i]; }                        //draw mixing prob values

    rdirichlet(prob,sumd,npat);

  }



}



void simnewsamples_ggC(double *xnew, int *dnew, double *anew, double *lnew, int *nsamples, int *groups, int *K, int *nsel, int *sel, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *gapprox) {

/* Simulates parameter values and new observations for a select subset of genes and the group indicated by groups from a Gamma/Gamma model

   Input arguments

   - nsamples: length of vector groups

   - groups: groups to simulate for.

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - nsel: expression values are simulated for nsel genes only.

   - sel: vector containing the keys of the nsel genes

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - balpha: estimate for b parameter in hyper-prior for alpha

   - v: matrix with post prob of each expression pattern (*nrow rows,*npat cols). In 2 group case reduces to post prob of equal and differential expression

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

  - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on.

  - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

  - nobsx: vector with number of terms in the sum. The first ngrouppat[0] elements correspond to pattern 0 and so on.

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used

   Output arguments

   - xnew: matrix of nsel rows and nsamples cols with obs drawn from the predictive for groups given by 'groups'.

   - dnew: matrix of nsel rows and nsamples cols with expression pattern indicators drawn from the posterior.

   - anew: matrix of nsel rows and nsamples cols with alpha parameters from the posterior for groups given by 'groups'.

   - lnew: matrix of nsel rows and nsamples cols with lambda parameters from the posterior for groups given by 'groups'.

*/



  int i, j, k, m, found, di, ci, ncolsumx, *colini, coldi, one=1, newton=2, usexpred=0, lengtha;

  double u, rsum, vcum, a, lambda, b1, b2, s, *vclus, *sumxpred, *prodxpred, *nobsxpred;



  sumxpred= dvector(0,1); prodxpred= dvector(0,1); nobsxpred= dvector(0,1); //never used: initialized to avoid warnings



  vclus= dvector(0,*nclust);

  colini= ivector(0,*npat);

  colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

  ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];



  for (i=0; i<(*nsel); i++) {                                 //generate expr. values for nsel genes only

    //printf("%d",i); //debug

    for (k=0; k<(*nsamples); k++) {                           //generate nsamples samples for each gene

      //draw expression pattern 

      di= 0; j=0; found=0; vcum= 0;                                              

      u= ranf();

      while ((found==0) && (j<(*npat-1))) {

        vcum += v[sel[i]*(*npat)+j];

        if (u<=vcum) { found= 1; di= j; }                  

        j++;

      }

      if (found==0) { di= (*npat) - 1; }                   //note: di=last pattern if it hasn't been assigned to another pattern yet



    // draw cluster indicator

      if (*nclust>1) {

	for (m=0, rsum=1; m<(*nclust); m++) {

	  vclus[m]= pdfcond_pat_clus(i,di,m,alpha0,nu,balpha,nualpha,ngrouppat,colini,ncolsumx,sumx,sumxpred,prodx,prodxpred,nobsx,nobsxpred,equalcv,usexpred) + log(rho[m]);

	  if (m>0) { vclus[m]= exp(vclus[m]-vclus[0]); rsum+= vclus[m]; }

	}

	vclus[0]= 1/rsum;

	for (m=1; m<(*nclust); m++) { vclus[m]= vclus[m]/rsum; }

	ci= rdisc(vclus,*nclust);

      } else {

	ci= 0;

      }



      //draw alpha,l 

      if (*equalcv == 1) {

	lengtha= ngrouppat[di];

	coldi= colini[di];

	b1= balpha[0]/nualpha[0] - prodx[sel[i]]; b2= alpha0[ci]/nu[ci];

	rcgammaC(&a,&one,nobsx+coldi,balpha,&b1,alpha0+ci,&b2,sumx+sel[i]*ncolsumx+coldi,&lengtha,&newton);

	lambda= 1.0/gengam(alpha0[ci]/nu[ci] + a*sumx[sel[i]*ncolsumx+coldi+patterns[di*(*K)+groups[k]]],alpha0[ci] + a*nobsx[coldi+patterns[di*(*K)+groups[k]]]);

      } else {

	lengtha= 1;

	coldi= colini[di]+patterns[di*(*K)+groups[k]];

	b1= balpha[0]/nualpha[0] - prodx[sel[i]*ncolsumx + coldi]; b2= alpha0[ci]/nu[ci];

	s= sumx[sel[i]*ncolsumx + coldi];

	rcgammaC(&a,&one,nobsx+coldi,balpha,&b1,alpha0+ci,&b2,&s,&lengtha,&newton);

	lambda= 1.0/gengam(alpha0[ci]/nu[ci] + a*s,alpha0[ci] + a*nobsx[coldi]);

      }



      //draw observation from the predictive

      xnew[i*(*nsamples)+k]= gengam(a/lambda,a);

      anew[i*(*nsamples)+k]= a;

      dnew[i*(*nsamples)+k]= di;

      lnew[i*(*nsamples)+k]= lambda;

    } //End k for

  }  //End i for



  free_ivector(colini,0,*npat);

  free_dvector(vclus,0,*nclust);

  free_dvector(sumxpred,0,1); free_dvector(prodxpred,0,1); free_dvector(nobsxpred,0,1);



}





void simpar_ggC(double *ngroupstot, double *sumd, double *sumci, double *sumalpha, double *sumlogalpha, double *suminvlambda, double *sumlambda, double *sumloglambda, int *groups, int *K, int *nrow, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *gapprox) {

/* Draws parameter values from the posterior, conditional on hyper-parameter values

   Input arguments

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - nrow: number of rows of x

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - balpha: estimate for b parameter in hyper-prior for alpha

   - nualpha: estimate for mu parameter in hyper-prior for alpha

   - v: matrix with post prob of each expression pattern (*nrow rows,*npat cols). In 2 group case reduces to post prob of equal and differential expression

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on. If usesumx==0 its values are computed but it still must be a vector of the appropriate size.

   - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

   - nobsx: vector with number of terms in the sum. The first ngrouppat[0] elements correspond to pattern 0 and so on.

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used



   Output arguments

   - ngroupstot: number of distinct (alpha,lambda) parameters

   - sumd: vector of length *npat indicating the number of genes falling into each expression pattern.

   - sumci: vector of length *nclust indicating the number of genes falling into each hyper-param cluster

   - sumalpha: sum of distinct alpha's (needed to draw from conditional posterior of hyper-parameters)

   - sumlogalpha: sum of distinct log(alpha)

   - suminvlambda: sum of 1/lambda

   - sumlambda: sum of distinct l's

   - sumloglambda: sum of distinct log l's

*/



  int i, j, m, found, di, ci=0, ncolsumx, *colini, one=1, newton=2, usexpred=0, lengtha;

  double u, vcum, rsum, lambda, a, b1, b2, s, *vclus, *sumxpred, *prodxpred, *nobsxpred;



  //FILE *ofile; //debug

  //ofile= openOut("simpar.txt"); //debug



  sumxpred= dvector(0,1); prodxpred= dvector(0,1); nobsxpred= dvector(0,1); //never used: initialized to avoid warnings



  vclus= dvector(0,*nclust);

  //Draw parameter values from the posterior

  for (i=0; i<(*npat); i++) { sumd[i]= 0; }

  for (m=0; m<(*nclust); m++) {

    ngroupstot[m]= 0;

    sumalpha[m]= sumlogalpha[m]= suminvlambda[m]= sumlambda[m]= sumloglambda[m]= 0; 

    sumci[m]= 0;

  }

  colini= ivector(0,*npat);

  colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

  ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];

  for (i=0; i<(*nrow); i++) {

    //draw expression pattern 

    di= 0; j=0; found=0; vcum= 0;                                              

    u= ranf();

    while ((found==0) && (j<(*npat-1))) {

      vcum += v[i*(*npat)+j];

      if (u<=vcum) { found= 1; di= j; }                  

      j++;

    }

    if (found==0) { di= (*npat) - 1; }                   //note: d[i]=last pattern if it hasn't been assigned to another pattern yet

    sumd[di] += 1;



    //fprintf(ofile,"%d ",di); //debug



    // draw cluster indicator

    if (*nclust>1) {

      for (m=0, rsum=1; m<(*nclust); m++) {

	vclus[m]= pdfcond_pat_clus(i,di,m,alpha0,nu,balpha,nualpha,ngrouppat,colini,ncolsumx,sumx,sumxpred,prodx,prodxpred,nobsx,nobsxpred,equalcv,usexpred) + log(rho[m]);

	if (m>0) { vclus[m]= exp(vclus[m]-vclus[0]); rsum+= vclus[m]; }

      }

      vclus[0]= 1/rsum;

      for (m=1; m<(*nclust); m++) { vclus[m]= vclus[m]/rsum; }

      ci= rdisc(vclus,*nclust);

    } else {

      ci= 0;

    }

    sumci[ci] += 1;

    ngroupstot[ci] += ngrouppat[di];  //number of distinct param values in cluster ci



    //fprintf(ofile,"%d ",ci); //debug



    // draw alpha and l for each group

    if (*equalcv == 1) { 

      lengtha= ngrouppat[di];

      b1= balpha[0]/nualpha[0] - prodx[i]; b2= alpha0[ci]/nu[ci];

      rcgammaC(&a,&one,nobsx+colini[di],balpha,&b1,alpha0+ci,&b2,sumx+i*ncolsumx+colini[di],&lengtha,&newton);

      sumalpha[ci] += a; sumlogalpha[ci] += log(a);

      //fprintf(ofile,"%f ", a); //debug

      for (j=0; j<ngrouppat[di]; j++) {

	lambda= 1.0/gengam(alpha0[ci]/nu[ci] + a*sumx[i*ncolsumx+colini[di]+j], alpha0[ci] + a*nobsx[colini[di]+j]);

	suminvlambda[ci] += 1.0/lambda; sumlambda[ci] += lambda; sumloglambda[ci] += log(lambda);

        //fprintf(ofile,"%f ", lambda); //debug

      }

      //if (di==0) fprintf(ofile,"%f \n",lambda); else fprintf(ofile,"\n"); //debug

    } else {

      lengtha= 1;

      for (j=0; j<ngrouppat[di]; j++) {

	b1= balpha[0]/nualpha[0] - prodx[i*ncolsumx + colini[di]+j]; b2= alpha0[ci]/nu[ci];

	s= sumx[i*ncolsumx + colini[di]+j];

	rcgammaC(&a,&one,nobsx+colini[di]+j,balpha,&b1,alpha0+ci,&b2,&s,&lengtha,&newton);

	lambda= 1.0/gengam(alpha0[ci]/nu[ci] + a*s, alpha0[ci] + a*nobsx[colini[di]+j]);

	sumalpha[ci] += a; sumlogalpha[ci] += log(a);

	suminvlambda[ci] += 1.0/lambda; sumlambda[ci] += lambda; sumloglambda[ci] += log(lambda);

      }

    }

  }  //End i for

  free_ivector(colini,0,*npat);

  free_dvector(vclus,0,*nclust);



  free_dvector(sumxpred,0,1); free_dvector(prodxpred,0,1); free_dvector(nobsxpred,0,1);



  //fclose(ofile); //debug



}



void simpred_ggC(double *xpred, int *d, double *alpha, double *l, int *usel, int *deltat, int *groups, int *K, int *nsel, int *sel, int *nrow, int *ncol, double *x, double *alpha0, double *nu, double *balpha, double *nualpha, int *equalcv, int *nclust, double *rho, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *prodx, double *nobsx, int *usesumx, int *gapprox) {

/* Draws parameter values from the posterior and *deltat observations for *K groups from the predictive distribution of a Gamma/Gamma model, conditional on hyper-parameter values

   Input arguments

   - usel: if usel==0 values of d,alpha and l are sampled from the posterior. If usel==1 xpred is sampled conditional on given l,alpha.

   - deltat: number of observations per group to draw from the predictive distribution

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - nsel: number of observations to be used for the classification

   - sel: vector of length nsel indicating the row numbers that the observations in xnew correspond to in v, sumx, and nobsx.

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: matrix with observations used to fit the model

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - balpha: estimate for b parameter in hyper-prior for alpha

   - v: matrix with post prob of each expression pattern (*nrow rows,*npat cols). In 2 group case reduces to post prob of equal and differential expression

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

   - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on. If usesumx==0 its values are computed but it still must be a vector of the appropriate size.

   - prodx: analogous to sumx but contains products in log scale i.e. sums of logarithms

   - nobsx: vector with number of terms in the sum. The first ngrouppat[0] elements correspond to pattern 0 and so on.

   - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, nobsx. For usesumx==0 they're computed from the data x.

   - gapprox: for gapprox==0 Gamma with matching mean & var is used, for gapprox==1 a faster version is used



   Output arguments

   - xpred: if genxpred==1 matrix of *nrow rows and K*deltat colums with obs drawn from the predictive (cols in group order). xpred not changed for genxpred=0

   - d: if usel==0 returns indicators for expression pattern drawn from the posterior distribution (in the two-group case typically d[i]==0 indicates equal expression and d[i]==1 differential expression). If usel==1 the provided d is used to generate xpred.

   - alpha: if usel==0 returns matrix of *nrow rows, *K columns with alpha draws from the posterior conditional on d. Column i corresponds to group i. If usel==1 the provided alpha is used to generate xpred. If equalcv==1 a vector of *nrow rows is returned instead of a matrix.

   - l: if usel==0 returns matrix of *nrow rows, *K columns with lambda drawn from the posterior conditional on d, alpha. Column i corresponds to group i. If usel==1 the provided l is used to generate xpred.

*/



  int i, j, m, found, group, di, ci, ncolpred, ncolsumx, *colini, one=1, newton=2, usexpred=0, lengtha;

  double u, rsum, vcum, *lambda, *a, b1, b2, s, *vclus, *sumxpred, *prodxpred, *nobsxpred;



  sumxpred= dvector(0,1); prodxpred= dvector(0,1); nobsxpred= dvector(0,1); //never used: initialized to avoid warnings



  ncolpred= (*K)*(*deltat);                                                              //number of columns of xpred



  //Draw parameter values from the posterior

  if ((*usel==0)) {

    vclus= dvector(0,*nclust);

    colini= ivector(0,*npat);

    colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

    ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];



    if ((*usesumx)==0) { compute_sumxC(sumx,prodx,nobsx,equalcv,nsel,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one); }

    for (i=0; i<(*nsel); i++) {

      //draw expression pattern 

      d[sel[i]]= 0; j=0; found=0; vcum= 0;                                              

      u= ranf();

      while ((found==0) && (j<(*npat-1))) {

        vcum += v[sel[i]*(*npat)+j];

        if (u<=vcum) { found= 1; d[sel[i]]= j; }                  

        j++;

      }

      if (found==0) { d[sel[i]]= (*npat) - 1; }                   //note: d[sel[i]]=last pattern if it hasn't been assigned to another pattern yet

      di= d[sel[i]];



      // draw cluster indicator

      if (*nclust>1) {

	for (m=0, rsum=1; m<(*nclust); m++) {

	  vclus[m]= pdfcond_pat_clus(i,di,m,alpha0,nu,balpha,nualpha,ngrouppat,colini,ncolsumx,sumx,sumxpred,prodx,prodxpred,nobsx,nobsxpred,equalcv,usexpred) + log(rho[m]);

	  if (m>0) { vclus[m]= exp(vclus[m]-vclus[0]); rsum+= vclus[m]; }

	}

	vclus[0]= 1/rsum;

	for (m=1; m<(*nclust); m++) { vclus[m]= vclus[m]/rsum; }

	ci= rdisc(vclus,*nclust);

      } else {

	ci= 0;

      }



      // draw alpha and l for each group

      if (*equalcv == 1) {

	lengtha= ngrouppat[di];

        a= dvector(0, 1); lambda= dvector(0, ngrouppat[di]);

	b1= balpha[0]/nualpha[0] - prodx[sel[i]]; b2= alpha0[ci]/nu[ci];

	rcgammaC(a,&one,nobsx+colini[di],balpha,&b1,alpha0+ci,&b2,sumx+sel[i]*ncolsumx+colini[di],&lengtha,&newton);

        alpha[sel[i]*(*K)] = a[0];

        for (j=0; j<ngrouppat[di]; j++) { 

          lambda[j]= 1.0/gengam(alpha0[ci]/nu[ci] + a[0]*sumx[sel[i]*ncolsumx+colini[di]+j], alpha0[ci] + a[0]*nobsx[colini[di]+j]);

	  group= patterns[di * (*K) + j];

	  l[sel[i]*(*K)+j] = lambda[group];

	}

	free_dvector(a, 0, 1); free_dvector(lambda, 0, ngrouppat[di]);

      } else {

	lengtha= 1;

        a= dvector(0, ngrouppat[di]); lambda= dvector(0, ngrouppat[di]);

        for (j=0; j<ngrouppat[di]; j++) { 

  	  b1= balpha[0]/nualpha[0] - prodx[sel[i]*ncolsumx + colini[di]+j]; b2= alpha0[ci]/nu[ci];

	  s= sumx[sel[i]*ncolsumx + colini[di]+j];

	  rcgammaC(a+j,&one,nobsx+colini[di]+j,balpha,&b1,alpha0+ci,&b2,&s,&lengtha,&newton);

          lambda[j]= 1.0/gengam(alpha0[ci]/nu[ci] + a[j]*s, alpha0[ci] + a[j]*nobsx[colini[di]+j]);

	}

	for (j=0; j<(*K); j++) {

	  group= patterns[di * (*K) + j];

	  alpha[sel[i]*(*K)+j] = a[group];

	  l[sel[i]*(*K)+j] = lambda[group];

	}

	free_dvector(a, 0, ngrouppat[di]); free_dvector(lambda, 0, ngrouppat[di]);

      }

    }  //End i for

    free_ivector(colini,0,*npat);

    free_dvector(vclus,0,*nclust);

  }  //End draw parameter values from the posterior



  //Draw observations from the predictive

  for (i=0; i<(*nsel); i++) {

    m=0;

    for (j=0; j<((*deltat)*(*K)); j++) {

      xpred[sel[i]*ncolpred+j]= gengam(alpha[sel[i]*(*K)+m]/l[sel[i]*(*K)+m],alpha[sel[i]*(*K)+m]);

      if (((j+1)%(*deltat)) == 0) m++;

    }

  }  //End i for



  free_dvector(sumxpred,0,1); free_dvector(prodxpred,0,1); free_dvector(nobsxpred,0,1);



}





void simpred_oldggC(double *xpred, int *d, double *l, int *usel, int *deltat, int *groups, int *K, int *nrow, int *ncol, double *x, double *alpha, double *alpha0, double *nu, double *v, int *npat, int *patterns, int *ngrouppat, double *sumx, double *nobsx, int *usesumx) {

/* Simulates *deltat observations for *K groups from the predictive distribution of a Gamma/Gamma model

   Input arguments

   - usel: if uselpred==0 values of d,l and xpred are sampled from the posterior. If uselpred==1 xpred is sampled conditional on given l.

   - deltat: number of observations per group to draw from the predictive distribution

   - groups: vector indicating what group each column in x corresponds to

   - K: number of groups e.g. group 0 for control, group 1 for type A cancer, group 2 for type B cancer

   - nrow: number of rows of x

   - ncol: number of columns of x

   - x: matrix with observations used to fit the model

   - alpha: estimate for alpha parameter in Gamma/Gamma model (shape parameter for observation component)

   - alpha0: estimate for alpha0 parameter in Gamma/Gamma model (shape parameter for hyper-prior)

   - nu: estimate for nu parameter in Gamma/Gamma model (location parameter for hyper-prior)

   - v: matrix with post prob of each expression pattern (*nrow rows,*npat cols). In 2 group case reduces to post prob of equal and differential expression

   - npat: number of patterns e.g. number of hypothesis regarding the mean expression levels.

   - patterns: matrix indicating which groups are put together under each pattern (*npat rows, K cols).

   - ngrouppat: vector of length *npat indicating the number of groups with different expression levels for each pattern

                (it can actually be computed from the patterns argument but it's more efficient to pre-compute and pass it as an argument)

  - sumx: sums of columns of x for each gene and each group in each expression pattern. The first ngrouppat[0] columns correspond to pattern 0, the next ngrouppat[1] to pattern 1 and so on. If usesumx==0 its values are computed but it still must be a vector of the appropriate size.

  - nobsx: vector with number of terms in the sum. The first ngrouppat[0] elements correspond to pattern 0 and so on.

  - usesumx: usesumx==1 indicates to use the sufficient statistics stored in sumx, nobsx. For usesumx==0 they're computed from the data x.



   Output arguments

   - xpred: matrix of *nrow rows and K*deltat colums with obs drawn from the predictive. Columns are in group order.

   - d: if usel==0 returns indicators for expression pattern drawn from the posterior distribution (in the two-group case typically d[i]==0 indicates equal expression and d[i]==1 differential expression). If usel==1 d is ignored.

   - l: if usel==1 returns matrix of *nrow rows, *K columns with lambda drawn from the posterior conditional on the patterns indicated in *d. Column i corresponds to group i. If usel==0 l is ignored.

*/



  int i, j, p, found, group, di, ncolpred, ncolsumx, *colini, one=1, *sel, equalcv=1;

  double u, vcum, *lambda, *prodx;



  ncolpred= (*K)*(*deltat);                                                              //number of columns of xpred



  colini= ivector(0,*npat);

  colini[0]= 0; for (i=1; i<(*npat); i++) { colini[i]= colini[i-1] + ngrouppat[i-1]; }   //column of sumx at which each pattern starts

  ncolsumx= colini[(*npat)-1] + ngrouppat[(*npat)-1];

  prodx= dvector(0,(*nrow));

  if ((*usesumx)==0) { 

    sel= ivector(0,*nrow);

    for (i=0; i<(*nrow); i++) { sel[i]= i; }

    compute_sumxC(sumx,prodx,nobsx,&equalcv,nrow,sel,&ncolsumx,ncol,x,groups,K,npat,patterns,ngrouppat,&one);

    free_ivector(sel,0,*nrow);

  }



  if ((*usel==0)) {                                           //If l values have to be generated

    for (i=0; i<(*nrow); i++) {

      //draw expression pattern 

      d[i]= 0; j=0; found=0; vcum= 0;                                              

      u= ranf();

      while ((found==0) && (j<(*npat-1))) {

        vcum += v[i*(*npat)+j];

        if (u<=vcum) { found= 1; d[i]= j; }                  

        j++;

      }

      if (found==0) { d[i]= (*npat) - 1; }                   //note: d[i]=last pattern if it hasn't been assigned to another pattern yet

      di= d[i];



      // draw l for each group

      lambda= dvector(0, ngrouppat[di]);                                            

      for (j=0; j<ngrouppat[di]; j++) { 

        lambda[j]= gengam((*nu) + sumx[i*ncolsumx+colini[di]+j], (*alpha) * nobsx[colini[di]+j] + (*alpha0));

      }

      for (j=0; j<(*K); j++) {

        group= patterns[di * (*K) + j];

        l[i*(*K)+j] = lambda[group];

      }

      free_dvector(lambda, 0, ngrouppat[di]);

    }  //End i for

  }  //End If l values have to be generated



  //draw observations from the predictive

  for (i=0; i<(*nrow); i++) {

    p=0;

    for (j=0; j<((*deltat)*(*K)); j++) {

      xpred[i*ncolpred+j]= gengam(l[i*(*K)+p],(*alpha));

      if (((j+1)%(*deltat)) == 0) p++;

    }

  }  //End i for



  free_ivector(colini,0,*npat);

  free_dvector(prodx,0,(*nrow));



}





/*********************************************************************************************

                            GAMMA SHAPE PARAMETER CONJUGATE PRIOR

*********************************************************************************************/



void dcgammaC(double *y, double *normk, double *x, int *n, double *a, double *b, double *c, double *d, double *r, double *s, int *lengtha, int *newton) {

  /* Returns density of a conjugate gamma shape distribution evaluated at x */

  /* Un-normalized density is (x^(p+1)/(r+s*x^(p+1)))^(a*x+d) y^(b-(p+1)*d-1) exp(-y*c) gamma(a*x+d) / (gamma(x)^a) */

  /* Input:

     - x: vector with points to evaluate the density at

     - n: length of x and y

     - a,b,c,d,r,s,p: parameters

     - lengtha: length of vectors a, s.

     Output: 

     - y: density of a conjugate gamma shape distribution with parameters a,b,c evaluated at x

     - normk: normalization constant.

   */

  int i, logscale=0; double aest, best;



gapprox_par(&aest,&best,a,b,c,d,r,s,lengtha,newton);   //Find parameters of approximating Gamma



//Find normalization constant & approximate Gamma density

*normk= kcgammaC(a,b,c,d,r,s,lengtha,newton,&logscale);

for (i=0; i<(*n); i++) { y[i]= exp(aest*log(best) - gamln(&aest) + (aest-1)*log(x[i])-best*x[i]); }



}



double kcgammaC(double *a, double *b, double *c, double *d, double *r, double *s, int *lengtha, int *newton, int *logscale) {

/* Finds normalization constant for a conjugate gamma shape distribution */

/* Un-normalized density is (x/(r+s*x))^(a*x+d) y^(b-d-1) exp(-y*c) gamma(a*x+d) / (gamma(x)^a) */

/* Input:

     - a,b,c,d,r,s: parameters

     - lengtha: length of vectors a, s.

     - logscale: logscale==0 returns result in original scale, logscale==1 in log scale

   Output: normalization constant 

   Details: density is approximated with Gamma matching the location of the mode and the value of the log-density

*/



 double aest, best, m, ans; 



gapprox_par(&aest,&best,a,b,c,d,r,s,lengtha,newton);   //Find parameters of approximating Gamma



//Find approximate normalization constant

if (aest>1) { m= (aest-1)/best; } else { m= aest/best; }

ans= logcgammaf(m,a,*b,*c,*d,*r,s,*lengtha) - aest*log(best)+gamln(&aest) - (aest-1)*log(m)+best*m;

if (*logscale == 0) return(exp(ans)); else return(ans);



}



void mcgammaC(double *normk, double *m, double *v, double *a, double *b, double *c, double *d, double *r, double *s,  int *lengtha, int *newton) {

/* Computes mean & variance for a conjugate gamma shape distribution */

/* Input:

     - a,b,c,d,r,s,p: parameters

     - lengtha: length of vectors a, s.

   Output: 

   - normk: normalization constant.

   - m: mean of a conjugate gamma shape distribution with parameters a,b,d

   - v: variance of a conjugate gamma shape distribution with parameters a,b,d

*/

  int logscale=0; double aest, best;



gapprox_par(&aest,&best,a,b,c,d,r,s,lengtha,newton);   //Find parameters of approximating Gamma



//Compute approx moments & norm constant

*normk= kcgammaC(a,b,c,d,r,s,lengtha,newton,&logscale);

*m= aest/best;

*v= (*m)/best;



}



void rcgammaC(double *x, int *n, double *a, double *b, double *c, double *d, double *r, double *s, int *lengtha, int *newton) {

/* Generates random draws from a conjugate gamma shape distribution by approximating it with a Gamma */

/* Input:

   - n: number of random draws to generate (length of x)

   - a,b,c,d,r,s,p: parameters

   - lengtha: length of vectors a, s.

   Output:

   - x: vector of length n with the random draws

*/

  int i;

  double aest, best;



gapprox_par(&aest,&best,a,b,c,d,r,s,lengtha,newton);   //Find parameters of approximating Gamma

for (i=0; i<(*n); i++) { x[i]= gengam(best,aest); }  //Generate random variates



}



void gapprox_par(double *aest, double *best, double *a, double *b, double *c, double *d, double *r, double *s, int *lengtha, int *newton) {

//Find parameters of approximating Gamma

/* Input: 

     - a,b,c,d,r,s,p parameters. 

     - lengtha: length of vectors a, s.

     - newton: number of Newton steps to take to find the mode of the distribution.



   Output: aest, best are the parameters of the approximating Gamma

*/

  int i; double m, fp, fpp, step;



//Gamma approx based on Stiling's formula

if (*lengtha == 1) {

  if (*a == 0) {    //for a==0 && d==0 distribution is exactly Gamma

    *aest = *b - *d;

    *best = *c;

  } else {    

    *aest = *b + .5*(*a) - .5;

    *best = *c + (*a)*log((*s)/(*a));

  }

} else {

  for (i=0,*aest= *b,*best= *c; i<(*lengtha); i++) {

    *aest += .5*a[i] - .5;

    *best += a[i]*log(s[i]/a[i]);

  }

}



if (*aest <1.0e-10) *aest= 1.0e-10;

if (*best <1.0e-10) *best= 1.0e-10;



//Newton step to better locate the maximum of the distrib

for (i=0, step=1, fpp=-1; (i< *newton) && (*aest>1) && (fabs(step)>.01) && (fpp<0); i++) {

  m= (*aest-1)/(*best);

  fp= logcgammafp(m,a,*b,*c,*d,*r,s,*lengtha); 

  fpp= logcgammafpp(m,a,*b,*c,*d,*r,s,*lengtha);

  step= fp/fpp;

  if ((fpp<0) && (logcgammaf(m,a,*b,*c,*d,*r,s,*lengtha)<logcgammaf(m-step,a,*b,*c,*d,*r,s,*lengtha))) { 

    *aest = 1 - (m-step)*(m-step)*fpp; 

    *best = -(m-step)*fpp;

  }

}



}







//void gapprox_par(double *aest, double *best, double *a, double *b, double *c, double *d, double *r, double *s, int *lengtha, int *newton) {

//Find parameters of approximating Gamma

/* Input: 

     - a,b,c,d,r,s,p parameters. For newton==1 a Newton step is taken to find the mode of the distribution.

     - lengtha: length of vectors a, s.



   Output: aest, best are the parameters of the approximating Gamma

*/

/*  int i; double m, mnew, fp, fpp, step;



//Gamma approx based on Stiling's formula

if (*lengtha == 1) {

  if (*a == 0) {    //for a==0 && d==0 distribution is exactly Gamma

    *aest = *b - *d;

    *best = *c;

  } else {    

    *aest = *b + .5*(*a) - .5;

    *best = *c + (*a)*log((*s)/(*a));

  }

} else {

  for (i=0,*aest= *b,*best= *c; i<(*lengtha); i++) {

    *aest += .5*a[i] - .5;

    *best += a[i]*log(s[i]/a[i]);

  }

}



//Newton step to better locate the maximum of the distrib

if ((*newton == 1) && (*aest > 1)) {

  m= (*aest-1)/(*best);

  fp= logcgammafp(m,a,*b,*c,*d,*r,s,*lengtha); 

  fpp= logcgammafpp(m,a,*b,*c,*d,*r,s,*lengtha);

  step= fp/fpp;

  if ((fabs(step) > .05) && (fpp<0) && (logcgammaf(m,a,*b,*c,*d,*r,s,*lengtha)<logcgammaf(m-step,a,*b,*c,*d,*r,s,*lengtha))) { 

    mnew= m- step;

    fp= logcgammafp(mnew,a,*b,*c,*d,*r,s,*lengtha);

    fpp= logcgammafpp(mnew,a,*b,*c,*d,*r,s,*lengtha);

    m= mnew;

    step= fp/fpp;

    if (fpp<0) { *aest = 1 - mnew*mnew*fpp; *best = -mnew*fpp; }

  }

}



}

*/



//Log-density of conjugate Gamma and its derivaties

double logcgammaf(double x, double *a, double b, double c, double d, double r, double *s, int lengtha) {

  int i; double t1, ans;

  for (i=0, ans=(b-d*(.0+lengtha)-1)*log(x) - x*c; i<lengtha; i++) {

    t1= a[i]*x+d;

    ans+= gamln(&t1)- a[i]*gamln(&x) + (a[i]*x+d)*(log(x/(r+s[i]*x)));

  }

  return(ans);

 }

double logcgammafp(double x, double *a, double b, double c, double d, double r, double *s, int lengtha) { 

  int i; double t1, ans;

  for (i=0, ans=(b-d*(.0+lengtha)-1)/x - c; i<lengtha; i++) {

    t1= a[i]*x+d;

    ans+= a[i]*(digamma(t1)-digamma(x)+log(x/(r+s[i]*x))) + t1*r/(r*x+s[i]*x*x);

  }

  return(ans);

}

double logcgammafpp(double x, double *a, double b, double c, double d, double r, double *s, int lengtha) { 

  int i; double t1, ans;

  for (i=0, ans= -(b-d*(.0+lengtha)-1)/(x*x); i<lengtha; i++) {

    t1= a[i]*x+d;

    ans+= a[i]*(a[i]*trigamma(t1)-trigamma(x))+r/(r*x+s[i]*x*x)*(2*a[i]-(a[i]*x+d)*(r+s[i]*2*x)/(r*x+s[i]*x*x));

  }

  return(ans);

}







