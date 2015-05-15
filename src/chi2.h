#ifndef __MYLIB_CHI2_H__
#define __MYLIB_CHI2_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct chi2_parameters {
  size_t n;
  double *x;
  double *y;
  double *sigma;
  double (*model_f) (double x, const gsl_vector *par);
  double (*model_df) (unsigned int i, double x, const gsl_vector *par);
  size_t npars;
} chi2_parameters;

int chi_f (const gsl_vector *X, void *p, gsl_vector *f);

int chi_df (const gsl_vector *X, void *par, gsl_matrix * J);

int chi_fdf (const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix * J);

#endif
