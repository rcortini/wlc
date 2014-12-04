#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "f_min.h"
#include "f_root.h"
#include "f_deriv.h"
#include "wlc.h"



/****************************************************************
 * EXACT FORMULAE
 ***************************************************************/



/* this function calculates the function to minimize to calculate the free energy
 * per unit length of a WLC. Equation 14 of Marko1995. */
double wlc_g_min_handle (double x, void *params) {
  double *p = (double *) params;
  double ex;
  double F = p [0];
  double lpb = p [1];
  gsl_complex z;
  GSL_SET_COMPLEX (&z, 4*x, 0);
  gsl_complex zex = gsl_complex_coth (z);
  ex = GSL_REAL (zex);
  return (x/(2*lpb) - F)*(ex - 1/(2*x));
}

/* calculates the free energy per unit length of a WLC, as a function of the applied end force */
double wlc_g_F (double F, double lpb) {
  int minimizer_result, iter, max_iter = 100;
  double x_lo, x_hi, x0;
  double fx_lo, fx_hi, fx0;
  double x_min, p [2];
  struct f_min_params fminp;

  /* we are not able to calculate values higher than this one */
  /* if (F>=WLC_F_MAX) return WLC_G_MAX; */

  /* sets all the parameters for the minimizer */
  x_lo = 1e-10;
  x_hi = 10.;
  x0 = (x_hi-x_lo)/2.;

  /* assigns the parameter to the function to minimize */
  p [0] = F;
  p [1] = lpb;

  /* see if the minimum and maximum chosen values are ok */
  fx_hi = wlc_g_min_handle (x_hi, p);
  fx0 = wlc_g_min_handle (x0, p);
  fx_lo = wlc_g_min_handle (x_lo, p);

  /* if the values are not okay, enter this cycle to fix them */
  iter = 0;
  while (fx0>fx_lo || fx0>fx_hi) {
    if (fx0>fx_hi)
      x_hi = x_lo + (x_hi-x_lo)*2.;
    else
      x_hi = (x_hi-x_lo)/2.;

    x0 = (x_hi-x_lo)/2.;
    fx0 = wlc_g_min_handle (x0, p);
    fx_hi = wlc_g_min_handle (x_hi, p);
    iter++;
    if (iter>max_iter) {
      fprintf (stderr, "wlc_g_F: max_iter hit! F = %f\n", F);
      exit (EXIT_FAILURE);
    }
  };

  /* set the parameters of the minimizer */
  fminp.verbose = 0;
  fminp.max_iter = 1000;
  fminp.func = wlc_g_min_handle;
  fminp.func_p = p;
  fminp.eps_abs = 1.e-4;
  fminp.eps_rel = 0.;
  fminp.type = gsl_min_fminimizer_brent;

  /* invoke the minimizer */
  minimizer_result = f_min (x_lo, x0, x_hi, &x_min, &fminp);

  /* check the result of the minimizer and return */
  if (minimizer_result==GSL_SUCCESS)
    return wlc_g_min_handle (x_min, p);
  else {
    fprintf (stderr, "wlc_g: minimization failed! F = %f\n", F);
    exit (EXIT_FAILURE);
  }
}

/* Helmholtz free energy of wlc as a function of rho */
double wlc_f_rho (double rho, double lpb) {
  double F = wlc_F_rho (rho, lpb);
  return wlc_g_F (F, lpb) + F*rho;
}

double wlc_g_F_handle (double F, void *p) {
  double *lpb = (double *) p;
  return wlc_g_F (F, *lpb);
}

/* rho (F) exact = d g_wlc (F)/dF */
double wlc_rho_F (double F, double lpb) {
  return -f_deriv (F, wlc_g_F_handle, &lpb);
}

double wlc_F_rho_handle (double F, void *params) {
  double *p = (double *) params;
  double lpb = p [0];
  double rho = p [1];
  return wlc_rho_F (F, lpb) - rho;
}

/* F (rho) exact is obtained inverting rho (F) exact */
double wlc_F_rho (double rho, double lpb) {
  int root_solver_ret_code;
  double p [2];
  double F, F_guess, x_lo, x_hi, fx_lo, fx_hi;
  struct f_root_params f_root_p;
  
  /* if force is too high */
  if (rho>=1.)
    return WLC_F_MAX;

  /* assigns parameters of the root solver */
  p [0] = lpb;
  p [1] = rho;
  f_root_p.f = wlc_F_rho_handle;
  f_root_p.p = p;
  f_root_p.verbose = 0;
  f_root_p.eps_rel = 0.;
  f_root_p.eps_abs = 1.e-5;
  f_root_p.max_iter = 100;
  f_root_p.type = gsl_root_fsolver_brent;

  /* bounding interval */
  F_guess = wlc_F_rho_interp (rho, lpb);
  /* TODO : correct these values */
  x_lo = F_guess/2.; /* these are values of FORCE */
  x_hi = F_guess*2.;

  /* check that we have a good initial interval */
  fx_lo = wlc_F_rho_handle (x_lo, p);
  fx_hi = wlc_F_rho_handle (x_hi, p);
  while (fx_lo*fx_hi>0.) {
    x_hi *= 2.;
    x_lo /= 2.;
    fx_lo = wlc_F_rho_handle (x_lo, p);
    fx_hi = wlc_F_rho_handle (x_hi, p);
  }

  /* calculate the root in the specified interval */
  root_solver_ret_code = f_root (x_lo, x_hi, &F, &f_root_p);

  /* check the status of the root solver and return */
  if (root_solver_ret_code==GSL_SUCCESS)
    return F;
  else {
    fprintf (stderr, "wlc_F_rho: root solver failed! rho = %f", rho);
    exit (EXIT_FAILURE);
  }
}


/****************************************************************
 * INTERPOLATION FORMULAE
 ***************************************************************/



/* this function expresses the relationship between force and relative extension
 * (rho = z/L) in the Worm-like Chain model in an interpolation. Formula due to 
 * J. Marko and E. Siggia, Macromolecules 1995 */
double wlc_rho_F_interp_handle (double rho, void *params) {
  double *p = (double * ) params;
  double F = p [0];
  double lpb = p [1];
  return rho + (1./((1.-rho)*(1.-rho)) - 1.)/4. - F*lpb;
}



/* this function calculates the value of rho, as a function of F */
double wlc_rho_F_interp (double F, double lpb) {
  int root_solver_ret_code;
  double p [2];
  double rho, x_lo, x_hi;
  struct f_root_params f_root_p;
  
  /* if force is zero, then rho is zero */
  if (F==0.)
    return 0.;

  /* assigns parameters of the root solver */
  p [0] = F;
  p [1] = lpb;
  f_root_p.f = wlc_rho_F_interp_handle;
  f_root_p.p = p;
  f_root_p.verbose = 0;
  f_root_p.eps_rel = 0.;
  f_root_p.eps_abs = 1.e-5;
  f_root_p.max_iter = 100;
  f_root_p.type = gsl_root_fsolver_brent;

  /* a quick glance returns the initial bounding interval */
  x_lo = 0.; /* these are values of RHO */
  x_hi = 1.-1.e-7;

  /* calculate the root in the specified interval */
  root_solver_ret_code = f_root (x_lo, x_hi, &rho, &f_root_p);

  /* check the status of the root solver and return */
  if (root_solver_ret_code==GSL_SUCCESS)
    return rho;
  else {
    fprintf (stderr, "wlc_rho_F: root solver failed! F = %f", F);
    exit (EXIT_FAILURE);
  }
}

/* returns the force of the WLC as a function of rho = z/L and the
 * value of the persistence length, in the Marko interpolation formula */
double wlc_F_rho_interp (double rho, double lpb) {
  if (rho>=1.)
    return WLC_F_MAX;
  else
    return (rho + (1./((1.-rho)*(1.-rho)) - 1.)/4.)/lpb;
}



/****************************************************************
 * HIGH FORCE LIMIT
 ***************************************************************/



/* free energy of WLC at high force */
double wlc_g_F_highforce (double F, double lpb) {
  return -F + sqrt (F/lpb);
}

/* rho (F) at high force */
double wlc_rho_F_highforce (double F, double lpb) {
  return 1.-1./(2.*lpb*sqrt (F/lpb));
}

/* F (rho) at high force */
double wlc_F_rho_highforce (double rho, double lpb) {
  double y = (1.-rho);
  return 1./(4.*lpb*y*y);
}

/* wlc_g at high force, as a function of rho */
double wlc_g_rho_highforce (double rho, void *p) {
  double *lpb = (double *) p;
  return -(2.*rho-1.)/(4*(*lpb)*(1.-rho)*(1.-rho));
}

/* derivative of the high force wlc_g wrt rho */
double wlc_g_rho_highforce_derivative (double rho, void *p) {
  double *par = (double *) p;
  double lpb = par [0];
  return -3*rho/(4*lpb*pow(1.-rho,3.));
}
