#pragma once
//! switch 2 adresses.
//! \param a
//! \param b
inline void switchAdr(double *&a,double *&b)
{
  double *c=a; a=b;b=c;
}
