
      Please, read the full documentation in Doc/ !

   We treat only homogeneous systems: dy/dt = f(y).


1. Adapt CMakeList to your needs. 

 For the Intel compiler: -DINTELCOMP : includes some pragmas.

2. Available examples are :

 - Oregonator.hpp (3 equations, stiff)

 - AvcF.hpp       (21 equations, stiff).

 - BZ.hpp         (a 1d PDE (the Belousov Zhabotinsky problem)
   Jacobian is sparse).
   
 - KPP.hpp        (a 1d PDE (the classical Fisher KPP problem).

   Modify "typedef" intructions in main.cc to choose one example.



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



 
