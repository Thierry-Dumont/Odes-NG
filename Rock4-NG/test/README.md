We solve a non linear problem, the Fisher-KPP equation in dimension 1
---------------------------------------------------------------------

du/dt= \nu u" + 0.5 u (1-u).

on the interval [0,1].

Thus the spectral radius of the linearized problem is bounded by 
4*\nu/h^2 + 1, where h=1/Size (Size is the number of unknowns, and we
discretize in space by unniform finite differences).



NB:
--

* add #define ROCK4_OMP in your main.cc if you want the code to run in
  parallel with openmp. You also need to parallelize your RHS (see
  FKPPOMP.h), and to
  define the number of threads you want to use by something like
  
  export OMP_NUM_THREADS= _number of threads_ .
  
* A file "result" is created, which contains the solution at time=xend.
  


