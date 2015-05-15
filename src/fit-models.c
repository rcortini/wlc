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
#include "fit-models.h"
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
  const size_t npars = x_init->size;
  nlin_fit_parameters fit_pars;
  multifit_results *fit_results = multifit_results_alloc (npars);

  /* initialize the fitter parameters */
  fit_pars.n = n;
  fit_pars.x = x;
  fit_pars.y = y;
  fit_pars.sigma = sigma;
  fit_pars.npars = npars;
  fit_pars.type = gsl_multifit_fdfsolver_lmsder;
  fit_pars.eps_abs = 1.e-4;
  fit_pars.eps_rel = 1.e-4;
  fit_pars.max_iter = 400;
  fit_pars.model_f = wlc_Marko_f;
  fit_pars.model_df = wlc_Marko_df;

  /* now fit */
  nlin_fit (x_init, &fit_pars, fit_results);

  /* print exit status */
  print_multifit_results (fit_results, 1);

  /* free memory and exit */
  multifit_results_free (fit_results);
  return 0;
}
