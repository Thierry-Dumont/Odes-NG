Odes-NG
=======

**C++ library of (Ordinary) Differential  Equations solvers.**

_(NG? because it is a  New Generation of solvers)._

This library contains:

1) optimized rewritings of classical **_stiff_** ODES solvers.

2) specialized versions of **_stabilized explicit solvers_**, mainly for  **_parabolic
PDEs_**.

_Stabilized explicit solvers (Rock*) are multi-threaded (openmp)._

It is yet experimental !
======================

Doxygen documentation:
---------------------

Just do;
 doxygen Doxyfile 

Documentation goes in Doc/

 * html/ the html documentation.

 * latex the full Latex documentation (type make to generate it).


The main page of the html documentation explains how to use these codes.


Structure:
---------

test/ contains  examples, and some documentation.

Subdirectories Fortran/ contains original programs from Hairer and Wanner.

sage/ subdirectory contains SageMath material used to build the codes.

Concept:
-------

User provides a _C++ class_ which should completely describe the problem:

* RHS (and if necessary Jacobian of RHS).
* Size (number of equations).
* Some operators if necessary.
* Constants.
* ...

This allow to inline functions in the case of a rather simple RHS.


