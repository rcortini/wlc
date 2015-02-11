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

#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "f_min.h"
#include "f_root.h"
#include "f_deriv.h"
#include "wlc.h"
#include "cavity.h"
#include "cavity_gradient.h"


/***************************************************************
 *                                                             *
 *  wlc library, by Ruggero Cortini (cortini@lptl.jussieu.fr)  *
 *                                                             *
 *  Formulas derived mainly from the research article          *
 *  J. Marko, E. Siggia, "Stretching DNA", Macromolecules 28   *
 *  (1995), 26: 8759--8770                                     *
 *  hereafter referred to as "Marko1995".                      *
 *                                                             *
 ***************************************************************/



/****************************************************************
 * EXACT FORMULAE
 ***************************************************************/



/* Function to minimize to obtain the Gibbs free energy per unit length.
 * Equation 14 of Marko1995. */
double wlc_g_min_handle (double x, void *params) {
  double *p = (double *) params;
  double ex;
  double F = p [0];
  double lpb = p [1];
  gsl_complex z, zex;
  GSL_SET_COMPLEX (&z, 4*x, 0);
  zex = gsl_complex_coth (z);
  ex = GSL_REAL (zex);
  return (x/(2*lpb) - F)*(ex - 1/(2*x));
}

/* Gibbs free energy per unit length, as a function of the applied end force */
double wlc_g_F (double F, double lpb) {
  int minimizer_result, iter, max_iter = 100;
  double x_lo, x_hi, x0;
  double fx_lo, fx_hi, fx0;
  double x_min, p [2];
  struct f_min_params fminp;

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
      wlc_error ("wlc_g_F: max_iter hit! F = %f\n", F);
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
    wlc_error ("wlc_g: minimization failed! F = %f\n", F);
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

/* rho (F) = -d g_wlc (F)/dF */
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
  int root_solver_ret_code, iter, max_iter = 100;
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
  x_lo = F_guess/2.; /* these are values of FORCE */
  x_hi = F_guess*2.;

  /* check that we have a good initial interval */
  fx_lo = wlc_F_rho_handle (x_lo, p);
  fx_hi = wlc_F_rho_handle (x_hi, p);
  iter = 0;
  while (fx_lo*fx_hi>0.) {
    x_hi *= 2.;
    x_lo /= 2.;
    fx_lo = wlc_F_rho_handle (x_lo, p);
    fx_hi = wlc_F_rho_handle (x_hi, p);
    iter++;
    if (iter>max_iter) {
      wlc_error ("wlc_F_rho: max_iter hit! rho = %f\n", rho);
      exit (EXIT_FAILURE);
    }
  }

  /* calculate the root in the specified interval */
  root_solver_ret_code = f_root (x_lo, x_hi, &F, &f_root_p);

  /* check the status of the root solver and return */
  if (root_solver_ret_code==GSL_SUCCESS)
    return F;
  else {
    wlc_error ("wlc_F_rho: root solver failed! rho = %f", rho);
    exit (EXIT_FAILURE);
  }
}

/* calculates the Gibbs free energy as a function of extension */
double wlc_g_rho (double rho, double lpb) {
  return wlc_g_F (wlc_F_rho (rho, lpb), lpb);
}



/****************************************************************
 * INTERPOLATION FORMULAE
 ***************************************************************/



/* interpolation formula: equation 7 of Marko1995 */
double wlc_rho_F_interp_handle (double rho, void *params) {
  double *p = (double * ) params;
  double F = p [0];
  double lpb = p [1];
  double y = (1.-rho);
  return rho + (1./(y*y) - 1.)/4. - F*lpb;
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
    wlc_error ("wlc_rho_F: root solver failed! F = %f", F);
    exit (EXIT_FAILURE);
  }
}

/* returns the force of the WLC as a function of rho = z/L and the
 * value of the persistence length, in the Marko interpolation formula */
double wlc_F_rho_interp (double rho, double lpb) {
  if (rho>=1.)
    return WLC_F_MAX;
  else {
    double y = 1.-rho;
    return (rho + (1./(y*y) - 1.)/4.)/lpb;
  }
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

/****************************************************************
 * CAVITY FRAMEWORK FOR POLYMER ELONGATION
 ***************************************************************/


/* compute cavity elongation rho as a function of force F */
double wlc_rho_F_cavity (double f, double bB, double JB){
  
    cavity_workspace *cav_w = (cavity_workspace *) cavity_workspace_alloc ();
    int i = 0;
    double *P, l=0., Z=0., integral;

    /* initialise cavity workspace */

    cavity_workspace_initialise (cav_w);
    
    /*iterate cavity equations to get the exact cavity marginal*/
    cavity_iterate_marginal_equations (cav_w, f, bB, JB);
    
    /* integrate z*t*P(t) over the unit sphere */
    for (P= cav_w -> marginal; P< cav_w -> marginal+__Ntheta__*__Nphi__; P++){
        
        /* evaluate the integral I(t) = exp(J * t*u) *P_c(u) */
        integral = cavity_integrate_marginal (cav_w, cav_w -> marginal, JB, i*__Ntheta__*__Nphi__);
        
        /* get the elongation = t*z * exp(b_B*f * t*z) * I(t)^2, Eq. (11) of Massucci et al. 2014 */
        l += *(cav_w -> cos_theta+ i/__Nphi__) * exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB) * integral*integral * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
        
        /* And increase the normalisation Z */
        Z += exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB) * integral*integral * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
        
        i++;
    }
    
    /* normalise the elongation */
    l /= Z;
    
    cavity_workspace_free (cav_w);
    return l;
}


/* compute cavity elongation rho as a function of force F and the gradient. */
/* This is also used to evaluate the susceptibility [Eq. (12), Massucci et al. (2014)] and the correlation length at fixed force [Eq. (13), Massucci et al. (2014)]*/
double wlc_rho_F_cavity_and_gradient (double f, double bB, double JB, double * drho_dbB, double * drho_dJB, double * xi_f){
    
    cavity_gradient_workspace *cav_w = (cavity_gradient_workspace *) cavity_gradient_workspace_alloc ();
    int i = 0;
    double *P, l=0., zeta_0=0., dl_dbB=0., dl_dJB=0., Z=0., dZ_dbB=0., dZ_dJB=0., integral, d_integral_dbB, d_integral_dJB, dZ, d_dZ_dbB, d_dZ_dJB, cos_theta;
    
    /* initialise the cavity workspace */

    cavity_gradient_workspace_initialise (cav_w);

    /*iterate cavity equations to get the exact cavity marginal*/
    cavity_iterate_gradient_marginal_equations (cav_w, f, bB, JB);
    
    /* integrate z*t*P(t) over the unit sphere */
    for (P= cav_w -> marginal; P< cav_w -> marginal+__Ntheta__*__Nphi__; P++){
        
        /* evaluate the integral I(t) = exp(J * t*u) *P_c(u) */
        integral = cavity_integrate_marginal_and_gradient_bB (cav_w, cav_w -> marginal, JB, i*__Ntheta__*__Nphi__);
        
        /* and its gradient wrt the parameters bB, JB */
        d_integral_dbB = cavity_integrate_marginal_and_gradient_bB (cav_w, cav_w -> d_marginal_dbB, JB, i*__Ntheta__*__Nphi__);
        
        d_integral_dJB = cavity_integrate_gradient_JB_marginal (cav_w, cav_w -> marginal, cav_w -> d_marginal_dJB, JB, i*__Ntheta__*__Nphi__);
        
        /* get the elongation = t*z * exp(b_B*f * t*z) * I(t)^2, Eq. (11) of Massucci et al. 2014 */
        
        dZ = exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB) * integral*integral * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__);
        
        d_dZ_dbB = exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB)*integral * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__) * ( *(cav_w -> cos_theta+ i/__Nphi__) * f * integral + 2. * d_integral_dbB) ;
        
        d_dZ_dJB = exp(f* *(cav_w -> cos_theta+ i/__Nphi__)*bB)*integral * *(cav_w -> w_cos_theta+ i/__Nphi__) * *(cav_w -> w_phi + i%__Nphi__) * 2. * d_integral_dJB;
        
        cos_theta = *(cav_w -> cos_theta+ i/__Nphi__);
        
        l +=  cos_theta * dZ;
        
        zeta_0 += cos_theta * cos_theta * dZ;
        
        dl_dbB += cos_theta * d_dZ_dbB;
        
        dl_dJB += cos_theta * d_dZ_dJB;
        
        
        /* And increase the normalisation Z */
        Z += dZ;
        dZ_dbB += d_dZ_dbB;
        dZ_dJB += d_dZ_dJB;
        
        i++;
    }
    
    /* normalise the elongation and properly calculate the connected self-correlation and the gradient */
    l /= Z;
    
    zeta_0 = zeta_0/Z - l * l;
    
    dl_dbB = dl_dbB/Z - l * dZ_dbB/Z;
    
    dl_dJB = dl_dJB/Z - l * dZ_dJB/Z;
    
    /* Assign the gradient to the workspace */
    
    cav_w -> dLdbB = dl_dbB;
    cav_w -> dLdJB = dl_dJB;
    
    /* And the correlation length at force f xi(f) [Eq. (13) of Massucci et al. 2014]*/
    
    cav_w -> xi = -bB/(log(1./f*dl_dbB - zeta_0) - log(1./f*dl_dbB + zeta_0));

    *drho_dbB = cav_w -> dLdbB;
    *drho_dJB = cav_w -> dLdJB;
    *xi_f = cav_w -> xi;

    cavity_gradient_workspace_free (cav_w);
    
    return l;
}
