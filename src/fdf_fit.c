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


/* print results of multidimensional fitting */
void print_multifit_results (multifit_results *fit_results, unsigned int vflag) {
  unsigned int i, j;
  const unsigned int dim = fit_results->dim;
  if (vflag) {
    printf ("MULTIDIMENSIONAL FITTING RESULTS:\n");
    printf ("\tGSL exit code = %s\n", gsl_strerror (fit_results->retcode));
    printf ("\tBest fit parameters:\n");
    for (i=0; i<dim; i++) {
      printf ("\tc [%d] = %.8e +/- %.8e\n", i, gsl_vector_get (fit_results->c, i), sqrt (gsl_matrix_get (fit_results->cov, i, i)));
    }
    printf ("\tCovariance matrix:\n");
    for (i=0; i<dim; i++) {
      printf ("\t ");
      for (j=0; j<dim; j++)
	printf ("%.8e ", gsl_matrix_get (fit_results->cov, i, j));
      printf ("\n");
    }
    printf ("\n\tchi2 = %g\n", fit_results->chisq);
  }
  else {
    for (i=0; i<dim; i++)
      printf ("%.8e %.8e\n", gsl_vector_get (fit_results->c, i), sqrt (gsl_matrix_get (fit_results->cov, i, i)));
  }
}



/* allocates memory for multidimensional fitting results */
multifit_results * multifit_results_alloc (unsigned int dim) {
  multifit_results *fit_results = (multifit_results *) malloc (sizeof (multifit_results));
  fit_results->cov = gsl_matrix_alloc (dim, dim);
  fit_results->c = gsl_vector_alloc (dim);
  fit_results->dim = dim;
  return fit_results;
}



/* frees the memory associated to the multidimensional fitting results */
void multifit_results_free (multifit_results *fit_results) {
  gsl_matrix_free (fit_results->cov);
  gsl_vector_free (fit_results->c);
  free (fit_results);
}



/* performs a non-linear best-fit of parameters from a set of weighted data and a model, and
 * its derivatives, based on the algorithm passed through the fit_p pointer, using a least-square minimization
 * method */
void nlin_fit (const gsl_vector *x_start, nlin_fit_parameters *fit_p, multifit_results *results) {
  int retcode;
  unsigned int iter = 0, npars = fit_p->npars, n = fit_p->n;
  const gsl_multifit_fdfsolver_type *T = fit_p->type;
  gsl_multifit_fdfsolver *s;
  gsl_multifit_function_fdf f;
  chi2_parameters chi2_p;

  /* init the chi2 parameters */
  chi2_p.model_f = fit_p->model_f;
  chi2_p.model_df = fit_p->model_df;
  chi2_p.n = n;
  chi2_p.npars = npars;
  chi2_p.x = fit_p->x;
  chi2_p.y = fit_p->y;
  chi2_p.sigma = fit_p->sigma;

  /* init function to fit: the chi2 of the user-specified function
   * versus the user-specified data. */
  f.f = &chi_f;
  f.df = &chi_df;
  f.fdf = &chi_fdf;
  f.n = n;
  f.p = fit_p->npars;
  f.params = &chi2_p;

  /* initialize the fitter */
  s = gsl_multifit_fdfsolver_alloc (T, n, npars);
  gsl_multifit_fdfsolver_set (s, &f, x_start);

  /* iterate */
  do {
    iter++;
    retcode = gsl_multifit_fdfsolver_iterate (s);
    if (retcode)
      break;
    retcode = gsl_multifit_test_delta (s->dx, s->x, fit_p->eps_abs, fit_p->eps_rel);
  }
  while (retcode == GSL_CONTINUE && iter < fit_p->max_iter);

  /* assign the fit vector and the covariance matrix */
  gsl_multifit_covar (s->J, 0.0, results->cov);
  gsl_vector_memcpy (results->c, s->x);

  /* the chi2 is the norm of the target function at the last
   * iteration */
  results->chisq = gsl_blas_dnrm2 (s->f);

  /* free memory and return */
  gsl_multifit_fdfsolver_free (s);
  results->retcode = retcode;
}



/* calculates the value of the chi^2, as a function of the best-fit vector of 
 * parameters, and the parameters of the function fitter */
double chi2_from_fit (gsl_vector *fit, nlin_fit_parameters *fit_pars) {
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
