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

#ifndef __FIT_H__
#define __FIT_H__

/* useful defines */
#define CHUNK_SIZE 10
#define MAX_LINE_SIZE 4096

#define FIT(i) gsl_vector_get(fit, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

#include <gsl/gsl_vector.h>

/* model wrappers */

/* Marko model */
double wlc_Marko_f (double z, const gsl_vector *par);
double wlc_Marko_df (unsigned int i, double z, const gsl_vector *par);

/* the fit function */
int wlc_Marko_fit (size_t n, double *x, double *y, double *sigma, gsl_vector *x_init);

#endif
