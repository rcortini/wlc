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

#ifndef __MYGSL_F_ROOT_H__
#define __MYGSL_F_ROOT_H__

#include <stdio.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

/* f_root.c */
struct f_root_params {
  unsigned int verbose;
  double (*f) (double x, void *p);
  void *p;
  const gsl_root_fsolver_type *type;
  double eps_abs;
  double eps_rel;
  unsigned int max_iter;
};

int f_root (double x_min, double x_max, double *root, struct f_root_params *par);

#endif
