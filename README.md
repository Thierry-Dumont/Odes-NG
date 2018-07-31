Odes-NG
=======


**C++ library of (Ordinary) Differential  Equations solvers, with routines adapted to Parabolic Equations.**

_(Why NG? because it is a  New Generation of solvers)._

This library contains:

1) optimized rewriting of classical **_stiff_** ODES solvers.

2) specialized versions of **_stabilized explicit solvers_**, mainly for  **_large systems resulting from parabolic
PDEs_**.

_Stabilized explicit solvers (Rock*) are multi-threaded (openmp)._

It is yet experimental !
======================

Numerical Methods:
-----------------

* Radau5-NG : rewriting of the classical Radau5 code (by Hairer and Wanner).

* Rock4-NG  : rewriting of the Rock4 code of A. Abdulle, with less memory copies.

* Rock4-L   : specialization of Rock4 for linear problems.

*  Rock2-NG: rewriting of the Rock2 code of A. Abdulle, with less
   memory copies. :new:


* Rodas-NG  : rewriting of the code Rodas (by Hairer and Wanner).

* SDIRKL    : SDIRK methods, for linear problems.

* SymplectikRK : Gaussian symplectic methods.

We hope to add a specialization of Rock2 for linear problems in the
near future.

Structure:
---------

* Odes-NG/common/include/: contains parts of code common to _all_ integrators.

* For each integrator, the directory contains:

   * include/ : the core integrator.
   * test/ : examples, and some documentation.
   * An optional sage/ sub-directory contains SageMath material used to build the codes.

Concept:
-------

User provides a _C++ class_ which should completely describe the problem:

* RHS (and if necessary Jacobian of RHS).
* Size (number of equations).
* Some operators if necessary.
* Constants.
* ...

This allow to inline functions in the case of a rather simple RHS.


The original codes for Radau5, Rodas, Rock2 and Rock4 can be found at
https://www.unige.ch/~hairer/software.html

Doxygen documentation:
---------------------

Just do:

```
 doxygen Doxyfile 
```
Documentation (html and latex) goes in Doc/

 The main page of the html documentation explains how to use these codes.

How to build and run the tests:
-------------------------------
```
cd Radau5-NG   (or Rock4-L, Rock4-NG, Rodas-NG and so on)
cd test
mkdir Build
cd Build
export CXX=your_compiler (with your compiler = g++ or icc or clang++; default:g++)
cmake ..
make
./run
```

Note that:
---------
* Odes-NG introduces the namespace **odes**.
* Rock* codes are multithreaded (OpenMP).

References:
----------
To learn about Ordinary Differential Equations solvers, you should read 
the _bible_:

* _Solving Ordinary Differential Equations I_, by Hairer,
Nørsett,, Wanner,
* _Solving Ordinary Differential Equations II Stiff and
  Differential-Algebraic Problems_ by Hairer and  Wanner,
* _Geometric Numerical Integration_ by Hairer, Lubich and  Wanner.

About Rock methods, you can also read
[this paper](https://mathicse.epfl.ch/files/content/sites/mathicse/files/Mathicse%20reports%202011/27.2011_AA.pdf).
