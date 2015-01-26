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
#include "cavity.h"
#include "cavity_alloc.h"
#include "cavity_init.h"
#include "cavity_scalar.h"
#include "cavity_integrate.h"

/******************************************************************
 *                                                                *
 *  cavity library, by Francesco A. Massucci                      *
 *          (francesco.massucci@urv.cat)                          *
 *                                                                *
 *  Formulas derived mainly from the research article             *
 *  F. A. Massucci, I. Perez Castillo, & C. J. Perez Vicente,     *
 *  "Cavity approach for modeling and fitting polymer stretching" *
 *  Physical Review E, 90, (2014), 5: 052708                      *
 *  hereafter referred to as "Massucci et al. (2014)".            *
 *                                                                *
 *****************************************************************/

/* integrate exp(J * t * u) * P_c (u) over u (i.e. angles theta and phi) */
/* integral in Eq. (9) of Massucci et al. 2014 */
double cavity_integrate_marginal (cavity_workspace *cav_w, double *p_c, double JB, int INIT) {
    
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
    
    
/* iterate Eq. (9) of Massucci et al. 2014 */
void cavity_iterate_marginal_equations (cavity_workspace *cav_w, double f, double bB, double JB){
    

    int i;
    double error, *P, *p1, *p2, Z, integral;
    
    /* The equations are iterated recursively. 2 arrays for the cavity marginal are used and swapped at each iteration */
    p1 = cav_w-> marginal;
    p2 = cav_w-> marginal_dummy;
    
        
    
    do {
        
        /* initialise the normalising factor and its derivative with respect to b_B and J_B */
        Z=0.;
        i=0;
        
        /* The cavity marginal depends on angles theta and phi */
        /* It is discretised in Ntheta x Nphi points */
            
        for(P = p2; P<p2+__Ntheta__*__Nphi__;P++){
            
            /* Perform the integrals of the cavity equations using Gaussian quadratures */
            integral = cavity_integrate_marginal (cav_w, p1, JB, i*__Ntheta__*__Nphi__);
            
            /* P = exp (b_B * f * z*t) * Integral(t) */
            *P = exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB)*integral;
            
            /* increase the normalization */
            Z += *P * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
            
            i++;
        }
        
        /* normalise the marginals and check for convergence */
        error = 0.;

        i = 0;
        
        for(P = p2; P<p2+__Ntheta__*__Nphi__;P++){
            
            *P /= Z;
            
            error+=fabs(*P-*(p1+i));
            
            i++;
            
        }
        
        /* swap p1, p2 for further iteration */
        P = p2;
        
        p2 = p1;
        
        p1 = P;
        
        }
        
        /* repeat until convergence */
        while(error>__TOL__);
}

/* compute cavity elongation rho as a function of force F */
double cavity_rho_F (double f, double bB, double JB){
  
    cavity_workspace *cav_w = (cavity_workspace *) cavity_workspace_alloc ();
    int i = 0;
    double *P, l=0., Z=0., integral;

    /* initialise cavity workspace */

    cavity_workspace_initialise (cav_w);
    
    /*iterate cavity equations to get the exact cavity marginal*/
    cavity_iterate_marginal_equations (cav_w, f, bB, JB);
    
    /* integrate z*t*P(t) over the unit sphere */
    for (P= cav_w -> marginal; P< cav_w -> marginal+__Ntheta__*__Nphi__; P++){
        
        /* evaluate the integral I(t) = exp(J * t*u) *P_c(u) */
        integral = cavity_integrate_marginal (cav_w, cav_w -> marginal, JB, i*__Ntheta__*__Nphi__);
        
        /* get the elongation = t*z * exp(b_B*f * t*z) * I(t)^2, Eq. (11) of Massucci et al. 2014 */
        l += *(cav_w -> cos_theta+ i/__Nphi__) * exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB) * integral*integral * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
        
        /* And increase the normalisation Z */
        Z += exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB) * integral*integral * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
        
        i++;
    }
    
    /* normalise the elongation */
    l /= Z;
    
    cavity_workspace_free (cav_w);
    return l;
}
