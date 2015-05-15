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

/* reads a vector of ncols columns of data (column number
 * specified in the cols vector) from a file */
unsigned int read_data (FILE *f_in, unsigned int ncols, unsigned int *cols, double ***data) {
  unsigned int i, j, n, vector_size;
  char word [MAX_LINE_SIZE];

  /* check the "cols" array */
  for (i=1; i<ncols; i++)
    if (cols[i]<=cols[i-1]) {
      wlc_error ("\"cols\" vector must be strictly increasing!\n");
      exit (EXIT_FAILURE);
    }

  /* initialize the vectors to read */
  vector_size = CHUNK_SIZE;
  *data = (double **) malloc (ncols * (sizeof (double *)));
  for (i=0; i<ncols; i++)
    (*data) [i] = (double *) malloc (vector_size * sizeof (double));

  /* scan the input file */
  n = 0;
  while (fgets (word, sizeof (word), f_in) != NULL) { 
    double val;

    /* expand the data array if necessary */
    if (n>vector_size-1) {
      vector_size += CHUNK_SIZE;
      for (i=0; i<ncols; i++)
	if (safe_realloc (vector_size, &(*data) [i])) {
	  wlc_error ("No more memory!\n");
	  exit (EXIT_FAILURE);
	}
    }

    /* if it is not a comment, scan a line */
    if (word [0] == '#')
      continue;
    else {
      int bytes_now=0, bytes_consumed=0;
      i=0;
      j=0;
      while (i<ncols) {
	if (sscanf (word+bytes_consumed, "%lf%n", &val, &bytes_now) == 1) {
	  bytes_consumed += bytes_now;
	  if (j==cols[i])
	    (*data) [i++] [n] = val;
	  j++;
	}
	else {
	  wlc_error ("Error reading column number %d!\n", j);
	  exit (EXIT_FAILURE);
	}
      }
      n++;
    }
  }

  /* return the number of lines read */
  return n;
}
