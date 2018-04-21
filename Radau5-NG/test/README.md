

   #### We treat only homogeneous systems: dy/dt = f(y).


1. Adapt CMakeList to your needs. 

 For the Intel compiler: -DINTELCOMP  includes some pragmas.

2. Available examples are :

 - Oregonator.hpp (3 equations, stiff)

 - AvcF.hpp       (21 equations, stiff).

 - BZ.hpp         (a 1d PDE (the Belousov Zhabotinsky problem)
   Jacobian is sparse).
   
 - KPP.hpp        (a 1d PDE (the classical Fisher KPP problem).

   Modify "typedef" instructions in main.cc to choose one example.

 We recommend to have a look at  _Oregonator.hpp_ which is quite simple.
In this example, the Jacobian is provided, but you can use it or not
by changing the 
flag _ComputeJacobianNumerically_ (if set to true, the  Jacobian is computed by
finite differences, and the method _Jacobian_ is not used).


3. define or not LOGRADAU5
    If defined, this will produce a "logfile" file.

4. Then:
```
mkdir Build
cd Build
cmake ..
make
./run 
```
to compile and execute the code.



 
