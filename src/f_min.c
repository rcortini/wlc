/* wlc, a simple library to calculate worm-like chain polymer functions
 *
 * Copyright (C) 2014, 2015  Ruggero Cortini

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

#include "f_min.h"

/* this function minimizes the one-dimensional function, starting at point
   x0, in the interval [x_min, x_max]. */
int f_min (double x_min, double x0, double x_max, double *minimum, struct f_min_params *par) {
  unsigned int iter;
  int status;
  double x_lo, x_hi;
  gsl_function func;

  /* allocates the memory for the minimizer */
  gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (par->type);

  /* defines the function to minimize and to be handed to the minimizer */
  func.function = par->func;
  func.params = par->func_p;

  /* initializes the minimizer */
  gsl_min_fminimizer_set (s, &func, x0, x_min, x_max);

  /* iterate */
  iter = 0;
  do {
    iter++;
    status = gsl_min_fminimizer_iterate (s);

    /* check interval for convergence */
    x_lo = gsl_min_fminimizer_x_lower (s);
    x_hi = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (x_lo, x_hi, par->eps_abs, par->eps_rel);

    if (par->verbose)
      printf ("iter %d: x_lo = %f x_hi = %f\n", iter, x_lo, x_hi);

  } while (status==GSL_CONTINUE && iter<par->max_iter);

  /* get the value of the minimum and free the memory */
  *minimum = gsl_min_fminimizer_x_minimum (s);
  gsl_min_fminimizer_free (s);

  return status;
}
