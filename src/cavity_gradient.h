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

#ifndef __CAVITY_GRADIENT_LIB_H__
#define __CAVITY_GRADIENT_LIB_H__

#include "cavity_macros.h"

typedef struct{
    
    /* arrays for 2 cavity marginals to be iterated*/
    double *marginal;
    double *marginal_dummy;
    
    double *d_marginal_dbB;
    double *d_marginal_dbB_dummy;
    
    double *d_marginal_dJB;
    double *d_marginal_dJB_dummy;
    
    /*arrays for cos(theta) and corresponding weight */
    double *cos_theta;
    double *w_cos_theta;
    
    /*arrays for phi and corresponding weight */
    double *phi;
    double *w_phi;
    
    /*arrays for scalar product values */
    double *scalar_prod;
    
    /* the persistence length at fixed force */
    double xi;
    
    /* gradient wrt the parameters bB, JB */
    double dLdbB;
    double dLdJB;
    
} cavity_gradient_workspace;

#include "cavity_gradient_alloc.h"
#include "cavity_gradient_init.h"

double cavity_integrate_marginal_and_gradient_bB (cavity_gradient_workspace *, double *, double, int);

double cavity_integrate_gradient_JB_marginal (cavity_gradient_workspace *, double *, double *, double, int);

void cavity_iterate_gradient_marginal_equations (cavity_gradient_workspace *, double, double, double);

#endif
