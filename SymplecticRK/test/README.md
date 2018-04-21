This program is extremely primitive and experimental!
====================================================
* In main.cc, choose your floatting point precision
```
  //typedef long double Double;
  typedef double Double;
```
Note that long double does not seem to be very safe, neither
normalized.

* Choose one of the problems:
```
  typedef Outer<Double> F;
  //typedef OuterUranus<Double> F;
  //typedef OuterNeptune<Double> F;
  //typedef Kepler2<Double> F;
```

Outer.hpp lets you simulate the solar system, where the inner planets
are accumulated with the sun.

* Then, compile and run thecode as explained in the main page of Odes-NG.
