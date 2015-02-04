/* utils.c
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

/* #include */
#include <stdio.h>
#include <stdarg.h>
#include "fit.h"
#include "wlc.h"

/* a function to display a message */
void wlc_message (char *text, ...) {
  va_list args;
  printf ("wlc: INFO: ");
  va_start (args, text);
  vfprintf (stdout, text, args);
  va_end (args);
  fflush (stdout);
}

/* a function to display an error message */
void wlc_error (char *text, ...) {
  va_list args;
  fprintf (stderr, "wlc: ERROR: ");
  va_start (args, text);
  vfprintf (stderr, text, args);
  va_end (args);
  fflush (stderr);
}

/* a function that checks the return value of fopen */
FILE *safe_fopen (const char *path, const char *mode) {
  FILE *fp = fopen (path, mode);
  if (fp == NULL) {
    wlc_error ("Could not open %s\n", path);
    exit (EXIT_FAILURE);
  }
  else
    return fp;
}

/* a function that checks the return value of realloc */
int safe_realloc (unsigned int new_vector_size, double **vector) {
  double *is_null;
  is_null = (double *) realloc (*vector, new_vector_size * sizeof (double));
  if (is_null == NULL)
    return 1;
  *vector = is_null;
  return 0;
}

/* reads data in <x> <y> <sigma> format from input file stream
 * and returns the number of data points that were read */
unsigned int read_data (const char *input_file, double **x, double **y, double **sigma) {
  unsigned int n, vector_size;
  char word [MAX_LINE_SIZE];
  FILE *f_in = safe_fopen (input_file, "r");

  /* scan the input file */
  vector_size = CHUNK_SIZE;
  n = 0;
  while (fgets (word, sizeof (word), f_in) != NULL) { 
    int n_vals;
    double z, F, stdv;

    /* expand the x, y and sigma arrays if necessary */
    if (n>vector_size-1) {
      vector_size += CHUNK_SIZE;
      
      if (safe_realloc (vector_size, x) ||
          safe_realloc (vector_size, y) ||
          safe_realloc (vector_size, sigma)) {
	wlc_error ("No more memory!\n");
	exit (EXIT_FAILURE);
      }
    }

    /* if it is not a comment, scan a line */
    if (word [0] == '#')
      continue;
    else {
      n_vals = sscanf (word, "%lf %lf %lf\n", &z, &F, &stdv);
      if (n_vals==3) {
	(*x) [n] = z;
	(*y) [n] = F;
	(*sigma) [n] = stdv;
      }
      else {
	wlc_error ("Incorrect input file format\n");
	exit (EXIT_FAILURE);
      }
      n++;
    }
  }
  fclose (f_in);

  return n;
}
