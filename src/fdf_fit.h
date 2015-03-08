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

#ifndef __MY_GSL_FDF_FIT__
#define __MY_GSL_FDF_FIT__

#include <gsl/gsl_multifit_nlin.h>

/* this structure contains all the information necessary to perform the
 * non-linear least-square fit */
struct fdf_fit_parameters {
  size_t n;
  double *x;
  double *y;
  double *sigma;
  double (*model_f) (double x, const gsl_vector *par);
  double (*model_df) (unsigned int i, double x, const gsl_vector *par);
  size_t p;
  struct data *d;
  double eps_abs;
  double eps_rel;
  size_t max_iter;
  const gsl_multifit_fdfsolver_type *type;
};

int chi_f (const gsl_vector *X, void *p, gsl_vector *f);

int chi_df (const gsl_vector *X, void *par, gsl_matrix * J);

int chi_fdf (const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix * J);

int fdf_fit (const gsl_vector *x_start, struct fdf_fit_parameters *fit_p, gsl_vector *fit, gsl_matrix *covar);

double chi2_from_fit (gsl_vector *fit, struct fdf_fit_parameters *fit_pars);

#endif
