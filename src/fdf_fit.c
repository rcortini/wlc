/* wlc, a simple library to calculate worm-like chain polymer functions
*
* Copyright (C) 2014, 2015  Ruggero Cortini, Francesco A. Massucci
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "fdf_fit.h"

/* this function calculates chi = F_i (x_i, y_i, p) = f(x_i, p) - y_i/sigma_i */
int chi_f (const gsl_vector *X, void *p, gsl_vector *f) {
  struct fdf_fit_parameters *fit_p = (struct fdf_fit_parameters *) p;
  size_t n = fit_p->n;
  double *x = fit_p->x;
  double *y = fit_p->y;
  double *sigma = fit_p->sigma;
  size_t i;

  for (i=0; i<n; i++) {
    double Yi = fit_p->model_f (x [i], X);
    gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
  }

  return GSL_SUCCESS;
}



/* calculates the derivative components of the chi function */
int chi_df (const gsl_vector *X, void *par, gsl_matrix * J) {
  struct fdf_fit_parameters *fit_p = (struct fdf_fit_parameters *) par;
  size_t n = fit_p->n, p = fit_p->p;
  double *sigma = fit_p->sigma;
  double *x = fit_p->x;
  size_t i, j;

  for (i=0; i<n; i++) {
    /* Jacobian matrix J(i,j) = dfi / dxj, */
    /* where fi = (Yi - yi)/sigma[i]       */
    double s = sigma [i];
    for (j=0; j<p; j++) {
      gsl_matrix_set (J, i, j, fit_p->model_df (j, x[i], X)/s); 
    }
  }
  return GSL_SUCCESS;
}



/* calculates chi function and its derivatives at the same time */
int chi_fdf (const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix * J) {
  chi_f (x, p, f);
  chi_df (x, p, J);

  return GSL_SUCCESS;
}



/* performs a non-linear best-fit of parameters from a set of weighted data and a model, and
 * its derivatives, based on the algorithm passed through the fit_p pointer, using a least-square minimization
 * method */
int fdf_fit (const gsl_vector *x_start, struct fdf_fit_parameters *fit_p, gsl_vector *fit, gsl_matrix *covar) {
  unsigned int iter = 0, p = fit_p->p, n = fit_p->n;
  const gsl_multifit_fdfsolver_type *T = fit_p->type;
  gsl_multifit_fdfsolver *s;
  int status;
  gsl_multifit_function_fdf f;

  /* init function to fit */
  f.f = &chi_f;
  f.df = &chi_df;
  f.fdf = &chi_fdf;
  f.n = n;
  f.p = fit_p->p;
  f.params = fit_p;

  /* initialize the fitter */
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, x_start);

  /* iterate */
  do
  {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s);

    /* printf ("status = %s\n", gsl_strerror (status)); */

    if (status)
      break;

    status = gsl_multifit_test_delta (s->dx, s->x, fit_p->eps_abs, fit_p->eps_rel);
  }
  while (status == GSL_CONTINUE && iter < fit_p->max_iter);

  /* assign the fit vector and the covariance matrix */
  gsl_multifit_covar (s->J, 0.0, covar);
  gsl_vector_memcpy (fit, s->x);

  /* printf ("status = %s\n", gsl_strerror (status)); */

  gsl_multifit_fdfsolver_free (s);

  return status;
}



/* calculates the value of the chi^2, as a function of the best-fit vector of 
 * parameters, and the parameters of the function fitter */
double chi2_from_fit (gsl_vector *fit, struct fdf_fit_parameters *fit_pars) {
  unsigned int i;
  double chi2 = 0.;
  gsl_vector *f = gsl_vector_alloc (fit_pars->n);
  chi_f (fit, fit_pars, f);

  for (i=0; i<fit_pars->n; i++) {
    double fi = gsl_vector_get (f, i);
    chi2 += fi*fi;
  }

  gsl_vector_free (f);
  return chi2;
}

