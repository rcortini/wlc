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

#ifndef __FDF_FIT_H__
#define __FDF_FIT_H__

#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include "chi2.h"

typedef struct multifit_results {
  gsl_vector *c;
  gsl_matrix *cov;
  unsigned int dim;
  int retcode;
  double chisq;
} multifit_results;

void print_multifit_results (multifit_results *fit_results, unsigned int vflag);

multifit_results * multifit_results_alloc (unsigned int degree);

void multifit_results_free (multifit_results *fit_results);

/* this structure contains all the information necessary to perform the
 * non-linear least-square fit */
typedef struct {
  size_t n;
  double *x;
  double *y;
  double *sigma;
  double (*model_f) (double x, const gsl_vector *par);
  double (*model_df) (unsigned int i, double x, const gsl_vector *par);
  size_t npars;
  double eps_abs;
  double eps_rel;
  size_t max_iter;
  const gsl_multifit_fdfsolver_type *type;
} nlin_fit_parameters;

void nlin_fit (const gsl_vector *x_start, nlin_fit_parameters *fit_p, multifit_results *results);

double chi2_from_fit (gsl_vector *fit, nlin_fit_parameters *fit_pars);

#endif
