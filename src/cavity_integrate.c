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

#include "cavity_integrate.h"

/* compute the weights for Gaussian-Legendre quadrature integration */
/* given integration endpoints x_low, x_high */
/* and number of points Npoints */
void abs_and_weights(double x_low, double x_high, double *x, double *w, int Npts){
    
    int half,i,j;
    
    double *abs, a, b, c, v, x_half, mid_range, factor;
    
    /* compute middle points */
    
    half = Npts/2 + (int) (Npts%2!=0);
    
    x_half = 0.5*(x_high+x_low);
    
    mid_range = 0.5*(x_high-x_low);
    
    i=0;
    
    /* compute abscissae and weights */
    for (abs = x; abs < x + half; abs++) {
        
        v=cos(M_PI*(i+0.75)/(Npts+0.5));
        
        do {
            
            a = 1.0;
            b = 0.0;
            
            for (j=0; j<Npts; j++) {
                
                c = b;
                
                b = a;
                
                a = ((2.0*j+1.0)*v*b-j*c)/(j+1);
                
            }
            
            v -= a*(v*v-1.0)/(Npts*(v*a-b));
            
        } while (fabs(a*(v*v-1.0)/(Npts*(v*a-b))) > __TOL__);
        
        /* assign values to abscissae and weights */
        *abs = x_half-mid_range*v;
        *(x +Npts-i-1) = x_half+mid_range*v;
        
        factor = (Npts*(v*a-b))/(v*v-1.0);
        
        *(w+ i) = 2.0*mid_range/((1.0-v*v)*factor*factor);
        *(w +Npts-i-1) = *(w+ i);
        
        i++;
    }
}
