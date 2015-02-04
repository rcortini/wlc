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

#ifndef __WLCLIB_H__
#define __WLCLIB_H__

/* Boltzmann constant in erg/K */
#define K_BOLTZMANN 1.380650400e-16

/* TODO: remove this */
#define WLC_F_MAX 255.900000
#define WLC_G_MAX -253.642700

#include <stdio.h>

/* exact formulae */
double wlc_g_F (double F, double lpb);

double wlc_f_rho (double rho, double lpb);

double wlc_rho_F (double F, double lpb);

double wlc_F_rho (double rho, double lpb);

double wlc_g_rho (double rho, double lpb);

/* interpolation formulae */
double wlc_g_F_interp (double F, double lpb);

double wlc_rho_F_interp (double F, double lpb);

double wlc_F_rho_interp (double rho, double lpb);

/* high force limit */
double wlc_g_F_highforce (double F, double lpb);

double wlc_g_rho_highforce (double rho, void *p);

double wlc_g_rho_highforce_derivative (double rho, void *p);

double wlc_rho_F_highforce (double F, double lpb);

double wlc_F_rho_highforce (double rho, double lpb);

/* utility functions defined in utils.c */
void wlc_message (char *text, ...);

void wlc_error (char *text, ...);

FILE *safe_fopen (const char *path, const char *mode);

int safe_realloc (unsigned int new_vector_size, double **vector);

unsigned int read_data (const char *input_file, double **x, double **y, double **sigma);

/* cavity routines for discrete models */
/* elongation with the cavity method */
double cavity_rho_F (double, double, double);

/* elongation gradient and persistence length with the cavity method */
double cavity_rho_F_and_gradient (double, double, double, double *, double *, double *);


#endif
