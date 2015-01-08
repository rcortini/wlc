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

#ifndef __MYGSL_F_MIN_H__
#define __MYGSL_F_MIN_H__

#include <stdio.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

/* f_min.c */
struct f_min_params {
  unsigned int verbose;
  double (*func) (double x, void *p);
  void *func_p;
  double eps_rel;
  double eps_abs;
  unsigned int max_iter;
  const gsl_min_fminimizer_type *type;
};

int f_min (double x_min, double x0, double x_max, double *minimum, struct f_min_params *par);

#endif
