 This implementation is for linear problems!
==========================================
		  
As example we solve the Heat equation, discretized in 1d with finite differences.

Available methods are:
---------------------

- RKMethod1 : Method described at page 100, table 6-5 in Hairer-Wanner T2, 2nd edition.
- RKMethod2 : Crouzeix method. cf. Hairer-Wanner page 100, figure 6.18.
- RKRS      : the Crank-Nicolson method.
- EulerImp  : the Implicit Euler method.

All these methods are considered as Diagonaly Implicit Runge-Kutta
methods.

To choose a method, just comment/uncomment the following lines in
main.cpp:

```
 //typedef RKMethod1 RK;
 //typedef RKMethod2 RK;
 //typedef RKRS RK;
 typedef CrN RK;
 //typedef EulerImp RK;
```

Note:
----
We use _SuperLU_ to solve sparse linear systems. You need
 to install it, or to replace it by something else. It is a standard
 package in popular Linux distributions.
