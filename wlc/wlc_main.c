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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "wlc.h"
#include "fit.h"

void print_usage (const char *program_name) {
  printf ("Usage: %s [-v] [-T <temperature>] <function> <function arguments>\n", program_name);
}

void print_help () {
  printf ("\n");
  printf ("wlc: a program to calculate worm-like chain semi-flexible polymer functions\n\n");
  printf ("Functions available:\n\n");
  printf ("Exact formulae:\n");
  printf ("\tF_rho <rho> <lpb>: the force as a function of relative extension\n");
  printf ("\trho_F <F> <lpb>: the relative extension as a function of force\n");
  printf ("\n");
  printf ("Cavity theory formulae:\n");
  printf ("\tcavity_rho_F <F> <bB> <JB>: the relative extension as a function of force\n");
  printf ("Options:\n");
  printf ("\t-v: verbose output\n");
  printf ("\t-h: print this help and exit\n");
  printf ("\t-T <temperature>: assign temperature in Kelvins\n");
  printf ("\t   (note: in this case all output will be given in pN or pN/nm)\n");
  printf ("Note:\n");
  printf ("\t -- values of persistence length should be given in nm\n");
  printf ("\t -- all output is given in units of kT/nm, unless -T option provided\n");
}

int main (int argc, char *argv []) {
  int c, vflag = 0, Tflag = 0;
  double T;
  char *function_name;
  const char *program_name = "wlc";

  if (argc < 2) {
    wlc_error ("Incorrect usage\n");
    print_usage (program_name);
    exit (EXIT_FAILURE);
  }

  /* parse command line options */
  while ((c = getopt (argc, argv, "T:v::h::")) != -1) {
    switch (c) {
      case 'v' :
	vflag = 1;
	break;
      case 'h' :
	print_usage (program_name);
	print_help ();
	exit (0);
	break;
      case 'T' :
	Tflag = 1;
	T = atof (optarg);
	break;
      default :
	print_usage (program_name);
	exit (EXIT_FAILURE);
    }
  };

  /* the first non-option argument must be the function name */
  function_name = argv [optind];

  /* now behavior depends on which function we wanted to calculate */
  if (strcmp (function_name, "F_rho")==0) {
    double rho, lpb, F;

    /* check that we have sufficient arguments */
    if (optind+2>=argc) {
      wlc_error ("Incorrect usage\n");
      print_usage (program_name);
      printf ("Usage: wlc [-v] F_rho <rho> <lpb>\n");
      exit (EXIT_FAILURE);
    }
    rho = atof (argv [optind+1]);
    lpb = atof (argv [optind+2]);
    F = wlc_F_rho (rho, lpb);

    /* if temperature was assigned, convert to pN */
    if (Tflag)
      F *= K_BOLTZMANN*T*1.e14;

    /* choose how output is given */
    if (vflag)
      printf ("rho = %.5e lpb = %.5e F = %.5e\n", rho, lpb, F);
    else
      printf ("%.5e\n", F);
  }
  else if (strcmp (function_name, "rho_F")==0) {
    double rho, lpb, F;

    /* check that we have sufficient arguments */
    if (optind+2>=argc) {
      wlc_error ("Incorrect usage\n");
      print_usage (program_name);
      printf ("Usage: wlc [-v] rho_F <F> <lpb>\n");
      exit (EXIT_FAILURE);
    }
    /* if temperature was assigned, convert to pN */
    F = atof (argv [optind+1]);
    lpb = atof (argv [optind+2]);

    if (Tflag)
      F /= (K_BOLTZMANN*T*1.e14);

    rho = wlc_rho_F (F, lpb);

    /* choose how output is given */
    if (vflag)
      printf ("F = %.5e lpb = %.5e rho = %.5e\n", F, lpb, rho);
    else
      printf ("%.5e\n", rho);
  }
  else if (strcmp (function_name, "cavity_rho_F")==0) {
    double rho, JB, bB, F;

    /* check that we have sufficient arguments */
    if (optind+3>=argc) {
      wlc_error ("Incorrect usage\n");
      print_usage (program_name);
      printf ("Usage: wlc [-v] cavity_rho_F <F> <bB> <JB>\n");
      exit (EXIT_FAILURE);
    }
    /* if temperature was assigned, convert to pN */
    F = atof (argv [optind+1]);
    bB = atof (argv [optind+2]);
    JB = atof (argv [optind+3]);

    if (Tflag)
      F /= (K_BOLTZMANN*T*1.e14);

    rho = cavity_rho_F (F, bB, JB);

    /* choose how output is given */
    if (vflag)
      printf ("F = %.5e bB = %.5e JB = %.5e rho = %.5e\n", F, bB, JB, rho);
    else
      printf ("%.5e\n", rho);
  }
  else if (strcmp (function_name, "cavity_rho_F_and_gradient")==0) {
    double rho, JB, bB, F, drho_dbB, drho_dJB, xi_f;

    /* check that we have sufficient arguments */
    if (optind+3>=argc) {
      wlc_error ("Incorrect usage\n");
      print_usage (program_name);
      printf ("Usage: wlc [-v] cavity_rho_F_and_gradient <F> <bB> <JB>\n");
      exit (EXIT_FAILURE);
    }
    /* if temperature was assigned, convert to pN */
    F = atof (argv [optind+1]);
    bB = atof (argv [optind+2]);
    JB = atof (argv [optind+3]);

    if (Tflag)
      F /= (K_BOLTZMANN*T*1.e14);

    rho = cavity_rho_F_and_gradient (F, bB, JB, &drho_dbB, &drho_dJB, &xi_f);

    /* choose how output is given */
    if (vflag)
      printf ("F = %.5e bB = %.5e JB = %.5e drho_dbB = %.5e drho_dJB = %.5e xi_f = %.5e rho = %.5e\n", F, bB, JB, drho_dbB, drho_dJB, xi_f, rho);
    else
      printf ("%.5e\n", rho);
  }
  else if (strcmp (function_name, "Marko_fit")==0) {
    int fit_result;
    unsigned int n;
    char *input_file;
    double lp0, L0;
    double *x, *y, *sigma;
    gsl_vector *x_init = gsl_vector_alloc (2);

    /* check that we have sufficient arguments */
    if (optind+3>=argc) {
      wlc_error ("Incorrect usage\n");
      print_usage (program_name);
      printf ("Usage: wlc Marko_fit <lp0> <L0> <input_file>\n");
      exit (EXIT_FAILURE);
    }

    /* get parameters */
    lp0 = atof (argv [optind+1]);
    L0 = atof (argv [optind+2]);
    input_file = argv [optind+3];

    /* alloc the first chunk of data */
    x = (double *) malloc (CHUNK_SIZE * sizeof (double));
    y = (double *) malloc (CHUNK_SIZE * sizeof (double));
    sigma = (double *) malloc (CHUNK_SIZE * sizeof (double));

    /* read data from input stream */
    n = read_data (input_file, x, y, sigma);

    /* fit data to chosen model */
    gsl_vector_set (x_init, 0, lp0);
    gsl_vector_set (x_init, 1, L0);
    fit_result = wlc_Marko_fit (n, x, y, sigma, x_init);

    /* TODO : put an error condition on fit_result */

    /* free memory */
    free (x);
    free (y);
    free (sigma);
    gsl_vector_free (x_init);

    return fit_result;
  }
  else {
    wlc_error ("Incorrect usage\n");
    print_usage (program_name);
    exit (EXIT_FAILURE);
  }

  return 0;
}
