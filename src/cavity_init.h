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

#ifndef __CAVITY_INIT_H__
#define __CAVITY_INIT_H__

#include <math.h>
#include "cavity_integrate.h"
#include "cavity_macros.h"
#include "cavity.h"
#include "cavity_scalar.h"

void cavity_initialise_marginal (cavity_workspace *);
void cavity_workspace_initialise (cavity_workspace *);

#endif
