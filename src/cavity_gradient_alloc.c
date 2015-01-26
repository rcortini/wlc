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

#include "cavity_gradient_alloc.h"

/* allocate the memory for the cavity workspace with its gradient arrays */
cavity_gradient_workspace *cavity_gradient_workspace_alloc (void){
    cavity_gradient_workspace *cav_wspace;
    
    cav_wspace = (cavity_gradient_workspace *) malloc (sizeof(cavity_gradient_workspace));
    
    /* allocate space for the cavity marginals and their derivatives to iterate */
    cav_wspace -> marginal = __ALLOC_MARGINAL__;
    cav_wspace -> marginal_dummy = __ALLOC_MARGINAL__;
    
    cav_wspace -> d_marginal_dbB = __ALLOC_MARGINAL__;
    cav_wspace -> d_marginal_dbB_dummy = __ALLOC_MARGINAL__;
    
    cav_wspace -> d_marginal_dJB = __ALLOC_MARGINAL__;
    cav_wspace -> d_marginal_dJB_dummy = __ALLOC_MARGINAL__;
    
    /* allocate cos(theta) and weights for the Gauss-Legendre integration */
    cav_wspace -> cos_theta = __ALLOC_COS_THETA__;
    
    cav_wspace -> w_cos_theta = __ALLOC_COS_THETA__;
    
    /* allocate phi and weights for the Gauss-Legendre integration */
    cav_wspace -> phi = __ALLOC_PHI__;
    cav_wspace -> w_phi = __ALLOC_PHI__;
    
    /* allocate memory for the scalar product */
    cav_wspace -> scalar_prod = __ALLOC_SCALAR_PRODUCT__;
    
    return cav_wspace;
}

/* free the memory of the cavity gradient workspace */
void cavity_gradient_workspace_free (cavity_gradient_workspace *cav_wspace){
    
    /* free the cavity marginals */
    free(cav_wspace -> marginal);
    free(cav_wspace -> marginal_dummy);

    free(cav_wspace -> d_marginal_dbB);
    free(cav_wspace -> d_marginal_dbB_dummy);
    
    free(cav_wspace -> d_marginal_dJB);
    free(cav_wspace -> d_marginal_dJB_dummy);
    
    
    /* free cos(theta) and weights */
    free(cav_wspace -> cos_theta);
    
    free(cav_wspace -> w_cos_theta);
    
    /* free phi and weights */
    free(cav_wspace -> phi);
    free(cav_wspace -> w_phi);
    
    /* free the scalar product */
    free(cav_wspace -> scalar_prod);
    
    /* free the workspace */
    free(cav_wspace);
}
