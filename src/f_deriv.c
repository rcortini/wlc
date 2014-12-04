#include "f_deriv.h"

/* this function calculates the numerical derivative of the function */
double f_deriv (double x, double (*func) (double x, void *p), void *func_p) {
  double dx, df, err;
  gsl_function f;
  double h = GSL_SQRT_DBL_EPSILON;

  /* set the fields for the function */
  f.function = func;
  f.params = func_p;

  /* set the increment of the function */
  dx = x==0. ? h : x*h;

  /* invokes the derivative calculator and returns */
  gsl_deriv_central (&f, x, dx, &df, &err);
  return df;
}
