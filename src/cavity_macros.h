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

#define EPS 1.e-5

#ifndef __Ntheta__

#define __Ntheta__ 15

#endif


#ifndef __Nphi__

#define __Nphi__ 15

#endif


#ifndef __TOL__

#define __TOL__ 1.e-5

#endif


#define __ALLOC_COS_THETA__ (double *) malloc(__Ntheta__*sizeof(double))

#define __ALLOC_PHI__ (double *) malloc(__Nphi__*sizeof(double))

#define __ALLOC_MARGINAL__ (double *) malloc(__Ntheta__*__Nphi__*sizeof(double))

#define __ALLOC_SCALAR_PRODUCT__ (double *) malloc(__Ntheta__*__Nphi__*__Ntheta__*__Nphi__*sizeof(double))
