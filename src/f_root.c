#include "f_root.h"

/* this function returns the solution to f (x) = 0 */
int f_root (double x_min, double x_max, double *root, struct f_root_params *par) {
  unsigned int iter;
  int status;
  double x_lo, x_hi;
  gsl_root_fsolver *s;

  /* defines the function to pass to the solver */
  gsl_function func;
  func.function = par->f;
  func.params = par->p;

  /* allocates and initializes the minimizer */
  s = gsl_root_fsolver_alloc (par->type);
  gsl_root_fsolver_set (s, &func, x_min, x_max);

  /* start the iteration to find the root */
  iter = 0;
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    *root = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, par->eps_abs, par->eps_rel);

    if (par->verbose) 
      fprintf (stderr, "%d: x = %f\n", iter, *root);
  }
  while (status==GSL_CONTINUE && iter<par->max_iter);

  /* free the memory and return the found root */
  gsl_root_fsolver_free (s);
  return status;
}
