#ifndef __MYGSL_F_ROOT_H__
#define __MYGSL_F_ROOT_H__

#include <stdio.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

/* f_root.c */
struct f_root_params {
  unsigned int verbose;
  double (*f) (double x, void *p);
  void *p;
  const gsl_root_fsolver_type *type;
  double eps_abs;
  double eps_rel;
  unsigned int max_iter;
};

int f_root (double x_min, double x_max, double *root, struct f_root_params *par);

#endif
