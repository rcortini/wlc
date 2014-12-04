/* wlc, a simple library to calculate worm-like chain polymer functions
 *
 * Copyright (C) 2014  Ruggero Cortini

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

#include "f_deriv.h"

/* this function calculates the numerical derivative of the function */
double f_deriv (double x, double (*func) (double x, void *p), void *func_p) {
  double dx, df, err;
  gsl_function f;
  double h = GSL_SQRT_DBL_EPSILON;

  /* set the fields for the function */
  f.function = func;
  f.params = func_p;

  /* set the increment of the function */
  dx = x==0. ? h : x*h;

  /* invokes the derivative calculator and returns */
  gsl_deriv_central (&f, x, dx, &df, &err);
  return df;
}
