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

#include "f_root.h"

/* this function returns the solution to f (x) = 0 */
int f_root (double x_min, double x_max, double *root, struct f_root_params *par) {
  unsigned int iter;
  int status;
  double x_lo, x_hi;
  gsl_root_fsolver *s;

  /* defines the function to pass to the solver */
  gsl_function func;
  func.function = par->f;
  func.params = par->p;

  /* allocates and initializes the minimizer */
  s = gsl_root_fsolver_alloc (par->type);
  gsl_root_fsolver_set (s, &func, x_min, x_max);

  /* start the iteration to find the root */
  iter = 0;
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    *root = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, par->eps_abs, par->eps_rel);

    if (par->verbose) 
      fprintf (stderr, "%d: x = %f\n", iter, *root);
  }
  while (status==GSL_CONTINUE && iter<par->max_iter);

  /* free the memory and return the found root */
  gsl_root_fsolver_free (s);
  return status;
}
