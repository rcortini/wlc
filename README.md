wlc
===

C library to compute worm-like chain polymer functions

Depends on the Gnu Scientific Library (GSL), available at
  http://www.gnu.org/software/gsl/

This library provides fast routines to calculate functions
related to the worm-like chain model of a semi-flexible polymer.

The functions are named
wlc_<A>_<x>_<regime>
where A (x) is calculated, in the regime <regime>.

Currently available:
  - g: Gibbs free energy as a function of force or relative extension
  - f: Helmholtz free energy as a function of relative extension
  - rho: relative extension as a function of force
  - F: force as a function of extension
  
Regimes currently available:
  - variational: using the variational formulae derived by Marko and Siggia
  - interpolation
  - high force

Formulas derived mainly from the research article
J. Marko, E. Siggia, "Stretching DNA", Macromolecules 28
(1995), 26: 8759--8770
hereafter referred to as "Marko1995".
 
