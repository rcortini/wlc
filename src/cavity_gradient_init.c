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

#include "cavity_gradient_init.h"

/* initialise (actual and dummy) cavity marginal */
void cavity_gradient_initialise_marginal (cavity_gradient_workspace *cav_w){
    
    double *p;
    
    for (p = cav_w -> marginal; p < cav_w -> marginal+__Ntheta__*__Nphi__; p++){
        
        /* initialise with uniform distribution in the unitary sphere */
        *p = 1./M_PI;
    }
    
    for (p = cav_w -> marginal_dummy; p < cav_w -> marginal_dummy+__Ntheta__*__Nphi__; p++){
        
        /* initialise with uniform distribution in the unitary sphere */
        *p = 1./M_PI;
    }
    
    for (p = cav_w -> d_marginal_dbB; p < cav_w -> d_marginal_dbB+__Ntheta__*__Nphi__; p++){
        
        /* initialise with uniform distribution in the unitary sphere */
        *p = 1./M_PI;
    }
    
    for (p = cav_w -> d_marginal_dbB_dummy; p < cav_w -> d_marginal_dbB_dummy+__Ntheta__*__Nphi__; p++){
        
        /* initialise with uniform distribution in the unitary sphere */
        *p = 1./M_PI;
    }
    
    for (p = cav_w -> d_marginal_dJB; p < cav_w -> d_marginal_dJB+__Ntheta__*__Nphi__; p++){
        
        /* initialise with uniform distribution in the unitary sphere */
        *p = 1./M_PI;
    }
    
    for (p = cav_w -> d_marginal_dJB_dummy; p < cav_w -> d_marginal_dJB_dummy+__Ntheta__*__Nphi__; p++){
        
        /* initialise with uniform distribution in the unitary sphere */
        *p = 1./M_PI;
    }
    
}


void cavity_gradient_workspace_initialise (cavity_gradient_workspace *cav_w){
    
    /* assign values and weights to be used in the Gauss-Legendre integration to cos(theta) and phi */
    abs_and_weights(-1,1, cav_w->cos_theta, cav_w->w_cos_theta,__Ntheta__);
    
    abs_and_weights(0,2*M_PI,cav_w->phi,cav_w->w_phi,__Nphi__);
    
    /* compute the scalar product t*u to be used within the integral in the different cavity routines */
    cavity_gradient_compute_scalar_products (cav_w);
    
    /* initialise the cavity marginal */
    cavity_gradient_initialise_marginal (cav_w);
}
