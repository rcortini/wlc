#ifndef __MYGSL_F_DERIV_H__
#define __MYGSL_F_DERIV_H__

#include <gsl/gsl_deriv.h>

double f_deriv (double x, double (*func) (double x, void *p), void *func_p);

#endif
