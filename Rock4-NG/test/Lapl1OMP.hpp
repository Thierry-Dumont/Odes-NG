#ifndef Lapl1OMP__h
#define Lapl1OMP__h
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include <iostream>
using namespace std;
class Lapl1
{
  const int Size;
  double nu,h,uh2,rspec;
  const double pi=4*atan(1.0);  
public:
  //! constructor
  //! \param _size  size of the system.
  Lapl1(int _size): Size(_size)
  {
    nu=0.01;
    h=1./(Size-1.e0); uh2=nu/(h*h); rspec=4.*uh2;
  }
  //! destructor.
  ~Lapl1(){} 
  //! y=F(x)
  //! \param x
  //! \param y
  inline void operator()(double * Restrict x,double * Restrict  y) const
  {
    ASSUME_ALIGNED(x); ASSUME_ALIGNED(y);
    y[0]= uh2*(x[1]-x[0]);
#include "Ivdep.hpp"
    for(int i=1;i<Size-1;i++)
       y[i]=uh2*(x[i-1]-2.0*x[i]+x[i+1]);
    y[Size-1]= uh2*(x[Size-2]-x[Size-1]);
  }
  //! return spectral radius.
  inline  double rho() const {
    return rspec;
  }
  //! initial conditions.
  //! \param x
  void init(double x[])
  {
    for(int i=0;i<Size;i++)
      x[i]=cos(4*pi*i*h);
  }
  //! return size.
  inline int size() const {return Size;}
};
#endif
