/* gamma lasso regression -- Matt Taddy 2013 */

#include <stdlib.h>
#include <string.h>
#include <Rmath.h>
#include <time.h>

#include "latools.h"
#include "rhelp.h"
#include "polysolve.h"
#include "lhdmod.h"

/* global variables */
int dirty = 0;

double *Y = NULL;
int n;
double ysum;
double *X = NULL;
int d;
int Nx;
int *xr = NULL;
int *xc = NULL;
int *xk = NULL;

double *E = NULL;
double *B = NULL;

int *mod = NULL; // 0 free, 1 log, 2 lasso, 3 ridge
double *par = NULL; // E[lambda], lambda, tau
double *rate = NULL;

double *H = NULL;
double *G = NULL;
double *D = NULL;

double *QN0 = NULL;
double *QN1 = NULL;
char *family = NULL;
double *tmpn;

double (*calcL)(int, double*, double*);
double (*calcG)(int, double*, int*, double*, double*);
double (*calcH)(int, double*, int*, double*, double);

/* global cleanup function */
void gamlr_cleanup(){
  if(!dirty) return;

  if(QN0){ free(QN0); QN0 = NULL; }
  if(QN1){ free(QN1); QN1 = NULL; }
  if(family){ free(family); family = NULL; }
  if(tmpn){ free(tmpn); tmpn = NULL; }

  dirty = 0;
}

/* negative log post */
double neglogpost(double *beta, double *eta)
{
  double L = calcL(n, eta, Y);
  int k;

  for(k=0; k<d; k++){
    if(mod[k]==-1)
      L += rate[k]*par[0]*log( rate[k] + fabs(beta[k]) );
	else if(mod[k]==1)
	  L += par[0]*rate[k]*(1.0 - log(par[0]*rate[k]/(rate[k]+fabs(beta[k]))));
    else if(mod[k]==2)
	  L += rate[k]*par[1]*fabs(beta[k]);
    else if(mod[k]==3)
      L += beta[k]*beta[k]*0.5*par[2];
  }
  return L;
}

/* updated beta and eta after parameter move */
void update(int k, double bchange)
{
  if(bchange == 0.0) return;
  for(int i=xk[k]; i<xk[k+1]; i++) E[xr[i]] += X[i]*(bchange);
  B[k] = B[k] + bchange;
}

/* the quadratic rooter */
double glroot(int k, double sgn)
{
  myassert(mod[k]==-1);

  double b, c, ghb;
  ghb = G[k]/H[k]-B[k];
  b = ghb + sgn*rate[k];
  c = rate[k]*(ghb*sgn + par[0]/H[k]);

  double *roots = solvequadratic(b,c);

  int ns = 0;
  double sln = 0.0;
  for(int h=1; h<=((int) roots[0]); h++)
    if(sgn == sign(roots[h])){ sln = roots[h]; ns++; }
  myassert(ns == 1);

  free(roots);
  return sln;
}

double imove(char fam)
{
  myassert(mod[0]==0);
  double dbet;
  switch( fam )
    {
    case 'l':
	  dbet = -G[0]/((double) n);
      break;
    case 'p':
	  dbet = log(ysum) - log(G[0]+ysum);
      break;
    case 'b':
	  dbet = -sign(G[0])*fmin(D[0],fabs(G[0])/H[0]);
      break;
    default: error( "unrecognized family type in imove." );
    }
  return dbet;
}

/* The gradient descent move for given direction */
double Bmove(int k, double sgn)
{
  double dbet,ghb,lambda;
  if(H[k] == 0.0) return -B[k]; // happens only if you have all zero predictors

  if(mod[k]==0){
	dbet = -G[k]/H[k]; // unpenalized
  }
  else if(mod[k]==-1){ // non-guaranteed gamma lasso
	ghb = G[k] - H[k]*B[k];
	if(fabs(ghb) < par[0]) dbet = -B[k];
	else dbet = glroot(k, -sign(ghb)) - B[k];
  }
  else if(mod[k]==1){ // guaranteed gamma lasso
	lambda = par[0]*rate[k]/(rate[k]+fabs(B[k]));
	ghb = G[k] - H[k]*B[k];
	if(fabs(ghb) < lambda) dbet = -B[k];
	else dbet = -(G[k]-sign(ghb)*lambda)/H[k];
  }
  else if(mod[k]==2){ // lasso
	lambda = rate[k]*par[1];
	ghb = G[k] - H[k]*B[k];
	if(fabs(ghb) < lambda) dbet = -B[k];
	else dbet = -(G[k]-sign(ghb)*lambda)/H[k];
  }
  else if(mod[k]==3){ // ridge
    dbet = -(G[k] + B[k]*par[2])/(H[k] + par[2]);
  }
  else error("unrecognized mod in Bmove.");

  if( (fabs(dbet) > D[k]) & (family[0] != 'l') ) dbet = sign(dbet)*D[k]; // trust region bounds
  return dbet;
}

/* Quasi-Newton update (meta: inlcudes previous qn steps)
   z0 is altered upon return as your qn solution  */
void QNmove(int dim, double *z0, double *z1, double *z2){
  double u,v,w;
  for(int j=0; j<dim; j++){
    u = z1[j]-z0[j];
    v = z2[j]-z1[j];
    if((u!=0) & (u!=v))
      { w = u/(u-v);
		z0[j] = (1.0-w)*z1[j] + w*z2[j]; }
    else z0[j] = z2[j];
  }
}

/*
 * Main Function: Rgamlr
 *
 * Coordinate descent for MAP estimation of penalized coefficients
 *
 */

void R_gamlr(int *famid, int *n_in, int *d_in, double *Y_in, double *ysum_in,
			 int *Nx_in, double *X_in, int *xr_in, int *xc_in, int *xk_in,
			 double *loadings, double *fitted,
			 int *mod_in, double *par_in, double *rate_in,
			 double *D_in,  double *G_in, double *H_in,
			 double *tol_in, double *quasinewton, int *verbalize)
{
  dirty = 1; // flag to say the function has been called
  time_t itime = time(NULL);  // time stamp for periodic R interaction
  family = (char*)  malloc(sizeof(char) * 9); // allocate family name string

  switch( *famid )
    {
    case 1:
	  strcpy(family, "linear");
      calcL = &lin_nllhd;
      calcG = &lin_grad;
      calcH = &lin_curvebound;
      break;
    case 2:
	  strcpy(family, "poisson");
	  calcL = &po_nllhd;
      calcG = &po_grad;
      calcH = &po_curvebound;
      break;
    case 3:
	  strcpy(family, "binomial");
      calcL = &bin_nllhd;
      calcG = &bin_grad;
      calcH = &bin_curvebound;
      break;
    default: error( "unrecognized family type." );
    }

  int i, k, t, verb;
  double tol, Bdiff, bchange;
  double qn, Lold, Lnew, Lqn, Ldiff;
  int dozero, numzero;

  /** Build everything **/
  verb = *verbalize;
  tol = *tol_in;
  n = *n_in;
  d = *d_in;

  tmpn = new_dzero(n);

  Y = Y_in;
  ysum = *ysum_in;

  X = X_in;
  Nx = *Nx_in;
  xr = xr_in;
  xc = xc_in;
  xk = xk_in;

  G = G_in;
  H = H_in;
  D = D_in;

  mod = mod_in;
  par = par_in;
  rate = rate_in;

  B = loadings;
  E = fitted;
  Lnew = neglogpost(B,E);
  if(isinf(Lnew) || isnan(Lnew))
    error(" Infinite or NaN initial fit ");

  qn = *quasinewton;
  if(qn>0.0) // quasi newton acceleration
    { QN0 = new_dvec(d); QN1 = new_dvec(d); }

  /* introductory print statements */
  if(verb)
    { myprintf(mystdout,
			   "*** %s regression for %d obsvervations with %d covariates ***\n",
			   family, n, d-1);
      myprintf(mystdout, "Objective L initialized at %g\n", Lnew); }

  /* optimize until objective stops improving */
  dozero = 1;
  Bdiff = Ldiff = tol*100.0;
  t = 0;
  while(Bdiff > tol & Ldiff > tol/1000.0){
    numzero = 0;
	Bdiff = 0.0;

    if(qn>0.0 & fabs(Bdiff) < qn*3.0)
      {	void *tmp = QN0; QN0 = QN1; QN1 = tmp; copy_dvec(QN1, B, d); }

    // loop through coefficients
    for(k=0; k<d; k++)
	  if( B[k] != 0.0 || dozero || mod[k]==0 || mod[k]==3 || t == 1 || (t+k)%10 == 0 ){

		// get current grad and curve
		G[k] = (*calcG)(xk[k+1]-xk[k], &X[xk[k]], &xr[xk[k]], E, Y);
		if(family[0]!='l' || t==0)
		  H[k] = (*calcH)(xk[k+1]-xk[k], &X[xk[k]], &xr[xk[k]], E, D[k]);

		// move
		if(k==0) bchange = imove(family[0]);
		else bchange = Bmove(k, sign(B[k]));

		if(bchange!=0.0){
		  // update
		  update(k,  bchange);
		  Bdiff += fabs(bchange);
		  // trust region
		  if(family[0] != 'l')
			D[k] = fmax(0.5*D[k], 2.0*fabs(bchange)); }

		// sum the zeros and check for escape from R
		if(B[k] == 0.0) numzero++;
		itime = my_r_process_events(itime); }

  // iterate
  t++;
  Lold = Lnew;
  Lnew = neglogpost(B,E);

  // accelerate
  if(qn>0 && fabs(Bdiff) < qn && t%3 == 0)
      {
		QNmove(d, QN0, QN1, B);
		zero_dvec(tmpn,n);
		for(i=0; i<Nx; i++) tmpn[xr[i]] += X[i]*QN0[xc[i]];
		Lqn = neglogpost(QN0, tmpn);
		if(Lqn < Lnew)
		  { if(verb) myprintf(mystdout, "QN log lhd jump: %g\n", Lnew-Lqn);
			Lnew = Lqn;
			copy_dvec(B, QN0, d);
			copy_dvec(E, tmpn, n); }
      }

    // nllhd diff
    Ldiff = Lold - Lnew;

	// print
    if(verb)
      myprintf(mystdout,
			   "t = %d: beta diff = %g (Ldiff = %g with %d zero loadings).\n",
			   Bdiff, t, Ldiff, numzero);

    // check for irregular exits
    if(Lnew!=Lnew || !isfinite(Lnew))
      error("non-finite likelihood. \n");
	if(Ldiff < 0.0 &&  fabs(Ldiff) > tol){
	  //warning("stopped due to non-monotonic convergence. \n");
	  Bdiff = 0; dozero = 0; }

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(Bdiff < tol)
      { Bdiff += tol+1.0;
		dozero = 1; }
  }

  /* normal exit */
  gamlr_cleanup();
}


void R_glselect(int *estlam, int *famid,
				int *n, double *Y, double *eta,
				int *p,  double *df,
				double *beta, double *hess, double *rate,
				double *mu,  double *LPY, double *LLHD)
{

  double nd = (double) *n;
  double pd = (double) *p;
  double m = *mu;

  /* likelihood */
  switch( *famid )
    {
    case 1:
	  LLHD[0] = -lin_nllhd(*n, eta, Y);
	  double sig2 = -2.0*LLHD[0]/(nd - df[0]);
	  LLHD[0] = LLHD[0]/sig2 - nd*log(sqrt(sig2));
	  m = m/sig2;
      break;
    case 2:
	  LLHD[0] = -po_nllhd(*n, eta, Y);
      break;
    case 3:
	  LLHD[0] = -bin_nllhd(*n, eta, Y);
      break;
    default: error( "unrecognized family type." );
    }

  /* null model */
  if((*p)==0)
	{ LPY[0] = LLHD[0]; return; }

  /* 2pi bit */
  double normconst = 0.5*pd*log(2.0*PI);

  /* log prior and log |Hessian| */
  double logdethess = 0.0;
  double logprior = pd*log(0.5*m);
  double r,s,b,h;
  for(int j=0; j<(*p); j++){
	r = rate[j]; s = m*r; b = fabs(beta[j]); h = hess[j];
	if(*estlam){
	  logprior += -(s+1)*log(1+b/r);
	  logdethess += log(h - s/((r+b)*(r+b)));
	}else{
	  logprior += log(r) - b*m*r;
	  logdethess += log(h); }
  }

  LPY[0] = normconst + LLHD[0] + logprior - 0.5*logdethess;

}
