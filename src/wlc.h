#ifndef __WLCLIB_H__
#define __WLCLIB_H__

/* TODO: remove this */
#define WLC_F_MAX 255.900000
#define WLC_G_MAX -253.642700

/* exact formulae */
double wlc_g_F (double F, double lpb);

double wlc_f_rho (double rho, double lpb);

double wlc_rho_F (double F, double lpb);

double wlc_F_rho (double rho, double lpb);

/* interpolation formulae */
double wlc_g_F_interp (double F, double lpb);

double wlc_rho_F_interp (double F, double lpb);

double wlc_F_rho_interp (double rho, double lpb);

/* high force limit */
double wlc_g_F_highforce (double F, double lpb);

double wlc_g_rho_highforce (double rho, void *p);

double wlc_g_rho_highforce_derivative (double rho, void *p);

double wlc_rho_F_highforce (double F, double lpb);

double wlc_F_rho_highforce (double rho, double lpb);

#endif
