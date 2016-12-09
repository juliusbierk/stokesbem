# Boundary Element Method for Stokes Flow
This is a very simple implementation (my first in fact) of the method of regularized boundary elements.
It uses very low order quadratures (6th), and does not add extra at the near-singular elements.
It has also not been thoroughly tested or used for any published work.

The script can use regularized Stokeslets or Blakelets. It accepts both surfaces and splines as inputs, see D. J. Smith *"A boundary element regularized Stokeslet method applied to cilia- and flagella-driven flow"*.

Requires C++ armadillo for solving.
