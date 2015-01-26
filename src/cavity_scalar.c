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

#include "cavity_scalar.h"

/* compute the scalar product t(theta, phi)*u(theta',phi') for all values of theta, phi, theta', phi' */
/* and store it into an array of size Ntheta^2*Nphi^2 for further use */
void cavity_compute_scalar_products (cavity_workspace *cav_w) {
    
    double *s, cos_theta1, cos_theta2, phi1, phi2, cos_phi1_phi2;
    int i=0;
    
    for (s= cav_w -> scalar_prod; s < cav_w -> scalar_prod+__Ntheta__*__Nphi__* __Ntheta__*__Nphi__; s++){
        
        /* compute cos(theta) and cos(theta') */
        cos_theta1 = *( cav_w -> cos_theta + i/(__Nphi__*__Ntheta__*__Nphi__));
        cos_theta2 = *( cav_w -> cos_theta + (i%(__Ntheta__*__Nphi__))/__Nphi__);
        
        /* compute phi and phi' */
        phi1 = *(cav_w -> phi + (i/(__Ntheta__*__Nphi__))%__Nphi__);
        phi2 = *(cav_w -> phi + (i%(__Ntheta__*__Nphi__))%__Nphi__);
        
        cos_phi1_phi2 = cos(phi1-phi2);
        
        *s = sqrt(1.-cos_theta1*cos_theta1)*sqrt(1.-cos_theta2*cos_theta2)*cos_phi1_phi2+cos_theta1*cos_theta2;
        
        i++;
        
    }
}
