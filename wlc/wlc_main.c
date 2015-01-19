/* wlc, a simple library to calculate worm-like chain polymer functions
 *
 * Copyright (C) 2014, 2015  Ruggero Cortini

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
  if (strcmp (function_name, "rho_F")==0) {
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
  else {
    wlc_error ("Incorrect usage\n");
    print_usage (program_name);
    exit (EXIT_FAILURE);
  }

  return 0;
}
