/* Model-specific likelihood functions */

#include <stdlib.h>
#include <Rmath.h>
#include "rhelp.h"
#include "lhdmod.h"

// Linear

double lin_nllhd(int n, double *e, double *y){
  double l = 0.0;
  for(int i=0; i<n; i++) l += 0.5*(y[i] - e[i])*(y[i] - e[i]);
  return l;
}

double lin_grad(int n, double *x, int *o, double *e, double *y){
  double g = 0.0;
  for(int i=0; i<n; i++) g += -x[i]*(y[o[i]]-e[o[i]]);
  return g;
}

double lin_curve(int n, double *x, int *o, double *e){
  // o,e are unused
  double h = 0.0;
  for(int i=0; i<n; i++) h += x[i]*x[i];
  return h;
}

double lin_curvebound(int n, double *x, int *o, double *e, double d){
  // o,e,d are unused
  return lin_curve(n, x, o, e);
}

// Poisson

double po_nllhd(int n, double *e, double *y){
  double l = 0.0;
  for(int i=0; i<n; i++) l += exp(e[i]) - y[i]*e[i];
  return l;
}

double po_grad(int n, double *x, int *o, double *e, double *y){
  double g = 0.0;
  for(int i=0; i<n; i++) g += -x[i]*(y[o[i]]-exp(e[o[i]]));
  return g;
}

double po_curve(int n, double *x, int *o, double *e){
  double h = 0.0;
  for(int i=0; i<n; i++) h += x[i]*x[i]*exp(e[o[i]]);
  return h;
}

double po_curvebound(int n, double *x, int *o, double *e, double d){
  double h = 0.0;
  for(int i=0; i<n; i++) h += x[i]*x[i]*exp(e[o[i]] + d*fabs(x[i]));
  return h;
}

// binomial

double bin_nllhd(int n, double *e, double *y){
  double l = 0.0;
  for(int i=0; i<n; i++) l += log(1 + exp(-y[i]*e[i]));
  return l;
}

double bin_grad(int n, double *x, int *o, double *e, double *y){
  double g = 0.0;
  for(int i=0; i<n; i++) g += -y[o[i]]*x[i]/(1.0 + exp(y[o[i]]*e[o[i]]) );
  return g;
}

double bin_curve(int n, double *x, int *o, double *e){
  double h = 0.0;
  for(int i=0; i<n; i++)
	h += x[i]*x[i]/(2.0 + exp(e[o[i]]) + 1.0/exp(e[o[i]]));
  return h;
}

double bin_curvebound(int n, double *x, int *o, double *e, double d){
  double h = 0.0;
  double E;
  for(int i=0; i<n; i++){
	E = exp(e[o[i]] - sign(e[o[i]])*d*fabs(x[i]));
	h += x[i]*x[i]/(2.0 + fmin(2.0, E + 1.0/E)); }
  return h;
}

// wrappers for R

void R_nllhd(int *modid, double *l, int *n, double *eta, double *Y){
  switch( *modid )
    {
    case 1:
      *l = lin_nllhd(*n, eta, Y);
	  break;
    case 2:
      *l = po_nllhd(*n, eta, Y);
	  break;
    case 3:
      *l = bin_nllhd(*n, eta, Y);
      break;
    default:
	  error( "unrecognized model type." );
    }
}

void R_grad(int *modid, double *G, int *d,
		   double *X, int *xr, int *xk,
		   double *eta, double *Y)
{
  double (*grad)(int, double*, int*, double*, double*);
  switch( *modid )
    {
    case 1:
      grad = &lin_grad; break;
    case 2:
      grad = &po_grad; break;
    case 3:
      grad = &bin_grad; break;
    default:
	  error( "unrecognized model type." );
    }
  for(int k=0; k<(*d); k++)
	G[k] = (*grad)(xk[k+1]-xk[k], &X[xk[k]], &xr[xk[k]], eta, Y);
}

void R_curve(int *modid, double *H, int *d,
			double *X, int *xr, int *xk,
			double *eta)
{
  double (*curve)(int, double*, int*, double*);
  switch( *modid )
    {
    case 1:
      curve = &lin_curve; break;
    case 2:
      curve = &po_curve; break;
    case 3:
      curve = &bin_curve; break;
    default:
	  error( "unrecognized model type." );
    }
  for(int k=0; k<(*d); k++)
	H[k] = (*curve)(xk[k+1]-xk[k], &X[xk[k]], &xr[xk[k]], eta);
}

