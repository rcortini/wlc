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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cavity_macros.h"
#include "cavity_gradient.h"
#include "cavity_gradient_alloc.h"
#include "cavity_gradient_init.h"
#include "cavity_gradient_scalar.h"
#include "cavity_integrate.h"

/******************************************************************
 *                                                                *
 *  cavity library with gradient routines, by                     *
 *  Francesco A. Massucci (francesco.massucci@urv.cat)            *
 *                                                                *
 *  Formulas derived mainly from the research article             *
 *  F. A. Massucci, I. Perez Castillo, & C. J. Perez Vicente,     *
 *  "Cavity approach for modeling and fitting polymer stretching" *
 *  Physical Review E, 90, (2014), 5: 052708                      *
 *  hereafter referred to as "Massucci et al. (2014)".            *
 *                                                                *
 *****************************************************************/

/* integrate exp(J * t * u) * f (u) over u (i.e. angles theta and phi) */
/* integral in Eq. (9) of Massucci et al. 2014 */
double cavity_integrate_marginal_and_gradient_bB (cavity_gradient_workspace *cav_w, double *p_c, double JB, int INIT) {
    
    double *s, integral=0.;
    int i=0;
    
    for (s=cav_w -> scalar_prod+INIT; s<cav_w -> scalar_prod+INIT+__Ntheta__*__Nphi__; s++){
        
        /* weights fort theta and phi are stored in w_theta and w_phi, respectively */
        /* i/__Nphi__ gives the (index of) angle theta and i%__Nphi__ the one of phi */
        integral += *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__) *exp(*s*JB) * *(p_c + i);
        
        i++;
        
    }
    
    return integral;
}

/* integrate exp(J * t * u) * df (u)/(d JB) over u (i.e. angles theta and phi) */
/* first integral in Eq. (15) of Massucci et al. 2014 */
double cavity_integrate_gradient_JB_marginal (cavity_gradient_workspace *cav_w, double *p_c, double *dp_c_db, double JB, int INIT) {
    
    double *s, integral=0.;
    int i=0;
    
    for (s=cav_w -> scalar_prod+INIT; s<cav_w -> scalar_prod+INIT+__Ntheta__*__Nphi__; s++){
        
        /* weights fort theta and phi are stored in w_theta and w_phi, respectively */
        /* i/__Nphi__ gives the (index of) angle theta and i%__Nphi__ the one of phi */
        integral += *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__) *exp(*s*JB) * (*s* (*(p_c + i)) + *(dp_c_db+i)) ;
        
        i++;
        
    }
    
    return integral;
}

/* Iterate the equations for the cavity marginals and their gradient wrt to the parameters bB, JB */
/* This solves Eqs. (9) and (15) of Massucci et al. (2014) */
void cavity_iterate_gradient_marginal_equations (cavity_gradient_workspace *cav_w, double f, double bB, double JB){
    
    
    int i;
    double error, *P, *p1, *p2, *dp1_dbB, *dp2_dbB, *dp1_dJB, *dp2_dJB, Z, dZ_dbB, dZ_dJB, integral, d_integral_dbB, d_integral_dJB;
    
    /* The equations are iterated recursively. 2 arrays for the cavity marginal are used and swapped at each iteration */
    /* the gradient is trated similarly */
    
    p1 = cav_w-> marginal;
    p2 = cav_w-> marginal_dummy;
    
    dp1_dbB = cav_w-> d_marginal_dbB;
    dp2_dbB = cav_w-> d_marginal_dbB_dummy;
    
    dp1_dJB = cav_w-> d_marginal_dJB;
    dp2_dJB = cav_w-> d_marginal_dJB_dummy;
    
    
    
    do {
        
        /* initialise the normalising factor and its derivative with respect to b_B and J_B */
        Z=0.;
        dZ_dbB =0.;
        dZ_dJB =0.;
        
        
        i=0;
        
        /* The cavity marginal depends on angles theta and phi */
        /* It is discretised in Ntheta x Nphi points */
        
        for(P = p2; P<p2+__Ntheta__*__Nphi__;P++){
            
            /* Perform the integrals of the cavity equations using Gaussian quadratures */
            integral = cavity_integrate_marginal_and_gradient_bB (cav_w, p1, JB, i*__Ntheta__*__Nphi__);
            d_integral_dbB = cavity_integrate_marginal_and_gradient_bB (cav_w, dp1_dbB, JB, i*__Ntheta__*__Nphi__);
            d_integral_dJB = cavity_integrate_gradient_JB_marginal (cav_w, p1, dp1_dJB, JB, i*__Ntheta__*__Nphi__);
            
            /* P = exp (b_B * f * z*t) * Integral(t) */
            *P = exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB)*integral;
            *(dp2_dbB+i) = exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB) * (f * integral* *(cav_w -> cos_theta+ i/__Nphi__) + d_integral_dbB);
            *(dp2_dJB+i) = exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB) * d_integral_dJB;
            
            /* increase the normalization */
            Z += *P * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
            dZ_dbB += *(dp2_dbB+i) * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
            dZ_dJB += *(dp2_dJB+i) * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
            
            i++;
        }
        
        /* normalise the marginals and check for convergence */
        /* of both the marginals and the gradient */
        error = 0.;
        
        i = 0;
        
        for(P = p2; P<p2+__Ntheta__*__Nphi__;P++){
            
            *P /= Z;
            
            *(dp2_dbB+i) = *(dp2_dbB+i)/Z - *P * dZ_dbB/Z;
            
            *(dp2_dJB+i) = *(dp2_dJB+i)/Z - *P * dZ_dJB/Z;
            
            error+=fabs(*P-*(p1+i)) + fabs(*(dp2_dbB+i)-*(dp1_dbB+i)) + fabs(*(dp2_dJB+i)-*(dp1_dJB+i));
            
            i++;
            
        }
        
        /* swap p1, p2 for further iteration */
        P = p2;
        
        p2 = p1;
        
        p1 = P;
        
        /* the same for the derivative wrt bB */
        P = dp2_dbB;
        
        dp2_dbB = dp1_dbB;
        
        dp1_dbB = P;
        
        /* and the derivative wrt JB */
        P = dp2_dJB;
        
        dp2_dJB = dp1_dJB;
        
        dp1_dJB = P;
        
        
    }
    
    /* repeat until convergence */
    while(error>__TOL__);
}
