#include "chi2.h"



/* this function calculates chi = F_i (x_i, y_i, p) = f(x_i, p) - y_i/sigma_i */
int chi_f (const gsl_vector *X, void *p, gsl_vector *f) {
  chi2_parameters *fit_p = (chi2_parameters *) p;
  size_t n = fit_p->n;
  double *x = fit_p->x;
  double *y = fit_p->y;
  double *sigma = fit_p->sigma;
  size_t i;

  for (i=0; i<n; i++) {
    const double Yi = fit_p->model_f (x [i], X);
    gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
  }

  return GSL_SUCCESS;
}



/* calculates the derivative components of the chi function */
int chi_df (const gsl_vector *X, void *par, gsl_matrix * J) {
  chi2_parameters *fit_p = (chi2_parameters *) par;
  size_t n = fit_p->n, npars = fit_p->npars;
  double *sigma = fit_p->sigma;
  double *x = fit_p->x;
  size_t i, j;

  for (i=0; i<n; i++) {
    /* Jacobian matrix J(i,j) = dfi / dxj, */
    /* where fi = (Yi - yi)/sigma[i]       */
    const double s = sigma [i];
    for (j=0; j<npars; j++) {
      gsl_matrix_set (J, i, j, fit_p->model_df (j, x[i], X)/s); 
    }
  }
  return GSL_SUCCESS;
}



/* calculates chi function and its derivatives at the same time, default method */
int chi_fdf (const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix * J) {
  chi_f (x, p, f);
  chi_df (x, p, J);
  return GSL_SUCCESS;
}
