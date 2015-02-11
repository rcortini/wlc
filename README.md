# wlc

C library to compute worm-like chain polymer functions, written by
- **Ruggero Cortini**, Université Pierre et Marie Curie, Paris, France
- **Francesco Alessandro Massucci**, Universitat Rovira i Virgili, Tarragona, Spain

This library provides fast routines to calculate functions
related to models of a semi-flexible polymer (generally, worm-like chain models).
Both a continuous model and a discrete model are considered.


## Dependencies

Depends on the Gnu Scientific Library (GSL), available at
  http://www.gnu.org/software/gsl/

## Usage

The functions are named
wlc_A_x_regime
where A (x) is calculated, in the regime "regime".

Currently available:
  - g: Gibbs free energy as a function of force or relative extension
  - f: Helmholtz free energy as a function of relative extension
  - rho: relative extension as a function of force
  - F: force as a function of extension
  
Regimes currently available:
  - variational: using the variational formulae derived by Marko and Siggia
  - interpolation
  - high force
  - cavity
  - cavity_gradient

Regime cavity and cavity_gradient solve the discrete model in F. A. Massucci et al (2014).
All other regimes are inherent to the model in J. Marko & E. Siggia (1995).

Provides also a program to quickly access to function values, named "wlc"

## References

Formulas for variational, interpolation, and high force regimes are taken from the research article:
- J. Marko, E. Siggia, "Stretching DNA", __Macromolecules__, 28, (1995), 26: 8759--8770
DOI: [10.1021/ma00130a008](http://dx.doi.org/10.1021/ma00130a008)

Cavity formulas are taken from the research article:
- F. A. Massucci, I. Perez Castillo, C. J. Perez Vicente, "Cavity approach for modeling and fitting polymer stretching", __Physical Review E__, 90, (2014), 5: 052708
DOI: [10.1103/PhysRevE.90.052708](http://dx.doi.org/10.1103/PhysRevE.90.052708)
 
