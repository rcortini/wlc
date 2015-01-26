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

#ifndef __CAVITY_LIB_H__
#define __CAVITY_LIB_H__

#include "cavity_macros.h"

/* a workspace structure to handle all memory needed by the cavity routines */
typedef struct{
    
    /* arrays for 2 cavity marginals to be iterated*/
    double *marginal;
    double *marginal_dummy;
    
    /*arrays for cos(theta) and corresponding weight */
    double *cos_theta;
    double *w_cos_theta;
    
    /*arrays for phi and corresponding weight */
    double *phi;
    double *w_phi;
    
    /*arrays for scalar product values */
    double *scalar_prod;
    
} cavity_workspace;

#include "cavity_alloc.h"
#include "cavity_init.h"

#endif
