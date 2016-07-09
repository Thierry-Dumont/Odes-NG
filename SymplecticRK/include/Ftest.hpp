#ifndef Ftest__h
#define Ftest__h
#include "OdesException.hpp"
namespace odes{
/////////////////////////////////////////////////////////////////////////
/// Perform stoping test vor different float representation.
/// We take a value (t>0) near the minimum representable number.
/// Iterations will be stoped if the tested value (=difference of 2
///  iterates) is lower than t, and 2 successive evaluations are no 
/// more decreasing.
////////////////////////////////////////////////////////////////////////
template<class Double> class Ftest
{
  // non specialized version. Call it => errr at compile time
};
/// the "double" (64 bits IEEE).
template<> class Ftest<double>
{
  double v,vold;
  double t;
public:
  Ftest()
  {
    v=0; vold=1.e+40;
    t=5.e-308;
    
  }
  inline long double value() const{return v;}
  inline void operator=(long double x)
  {
    vold=v; v=x;
  }
  inline void operator+=(long double x){v+=x;}
  inline bool success() const {
    return v<=t && v>=vold;
  }
};
// the long double version.
template<> class Ftest<long double>
{
  long double v,vold;
  long double t;
public:
  Ftest()
  {
    v=0.0; vold=1.e+40;
    t=1.e-310;
  }
  inline long double value() const{return v;}
  inline void operator=(long double x)
  {
    vold=v; v=x;
  }
  inline void operator+=(long double x){v+=x;}
  inline bool success() const {
    return v<t && v>=vold;
  }

};
};
#endif
