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

#include "wlc.h"
#include "fit.h"
#include "fdf_fit.h"

/* model wrappers */

/* Marko model */
double wlc_Marko_f (double z, const gsl_vector *par) {
  double lpb = gsl_vector_get (par, 0);
  double L = gsl_vector_get (par, 1);

  return wlc_F_rho (z/L, lpb);
}

/* derivative of Marko model */
double wlc_Marko_df (unsigned int i, double z, const gsl_vector *par) {
  double lpb = gsl_vector_get (par, 0);
  double L = gsl_vector_get (par, 1);

  if (i==0) {
    double rho = z/L;
    return rho + (1./((1.-rho)*(1.-rho)) - 1.)/4.;
  }
  else if (i==1) {
    double rho = z/L;
    return -rho/(L*lpb)*(1 + 0.5/((1.-rho)*(1.-rho)*(1.-rho)));
  }
  else {
    wlc_error ("Invalid i = %d\n", i);
    exit (EXIT_FAILURE);
  }
}

/* fits data to Marko model */
int wlc_Marko_fit (size_t n, double *x, double *y, double *sigma, gsl_vector *x_init) {
  int fit_result, exit_code;
  const size_t p = x_init->size;
  struct fdf_fit_parameters fit_pars;
  gsl_vector *fit = gsl_vector_alloc (p);
  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  /* initialize the fitter parameters */
  fit_pars.n = n;
  fit_pars.x = x;
  fit_pars.y = y;
  fit_pars.sigma = sigma;
  fit_pars.p = p;
  fit_pars.type = gsl_multifit_fdfsolver_lmsder;
  fit_pars.eps_abs = 1.e-4;
  fit_pars.eps_rel = 1.e-4;
  fit_pars.max_iter = 400;
  fit_pars.model_f = wlc_Marko_f;
  fit_pars.model_df = wlc_Marko_df;

  /* now fit */
  fit_result = fdf_fit (x_init, &fit_pars, fit, covar);

  /* print exit status */
  wlc_message ("fit status = %s\n", gsl_strerror (fit_result));

  /* print fit result if success */
  if (fit_result==GSL_SUCCESS || fit_result == GSL_CONTINUE) {
    double chi2 = chi2_from_fit (fit, &fit_pars);
    double dof = n-p;
    double c = GSL_MAX_DBL(1, sqrt(chi2/dof));
    wlc_message ("chisq/dof = %g\n",  chi2/dof);
    wlc_message ("lpb     = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    wlc_message ("L       = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    exit_code = 0;
  }
  else
    exit_code = 1;

  /* free memory and exit */
  gsl_vector_free (fit);
  gsl_matrix_free (covar);
  return exit_code;
}
