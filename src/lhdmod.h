#ifndef __LHDMOD_H__
#define __LHDMOD_H__

double lin_nllhd(int n, double *e, double *y);
double lin_grad(int n, double *x, int *o, double *e, double *y);
double lin_curve(int n, double *x, int *o, double *e);
double lin_curvebound(int n, double *x, int *o, double *e, double d);

double po_nllhd(int n, double *e, double *y);
double po_grad(int n, double *x, int *o, double *e, double *y);
double po_curve(int n, double *x, int *o, double *e);
double po_curvebound(int n, double *x, int *o, double *e, double d);

double bin_nllhd(int n, double *e, double *y);
double bin_grad(int n, double *x, int *o, double *e, double *y);
double bin_curve(int n, double *x, int *o, double *e);
double bin_curvebound(int n, double *x, int *o, double *e, double d);

void Rnllhd(int *modid, double *l, int *n, double *eta, double *Y);
void Rgrad(int *modid, double *G, int *d, double *X, int *xr, int *xk, double *eta, double *Y);
void Rcurve(int *modid, double *H, int *d, double *X, int *xr, int *xk, double *eta);

#endif
