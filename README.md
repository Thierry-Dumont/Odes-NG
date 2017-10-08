odes
====

**C++ library of (Ordinary) Differential  Equations solvers.**


This library contains:

1) optimized rewritings of classical **_stiff_** ODES solvers.

2) specialized versions of stabilized explicit solvers, mainly for parabolic
PDEs.


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

Fortran/ subdirectory contains original programs from Hairer and Wanner.

sage/ subdirectory contains SageMath material used to build the codes.