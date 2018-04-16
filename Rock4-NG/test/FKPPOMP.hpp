#ifndef FKPPOMP__h
#define FKPPOMP__h
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include <iostream>
using namespace std;
class FKPP
{
  // This is the Fisher-KPP equation in dimension 1.
  // du/dt=  \nu u" + 0.5 u (1-u).
  // 
  const int Size;
  double nu,h,uh2,rspec;
  double f(double u) const 
  {
    // bound f in [0,1] is necessary due to floating approximations 
    // (otherwise, we will surely return values >1 => explosion !).
    return min(1.,max(0.,0.5*u*(1-u)));

  }
public:
  //! constructor
  //! \param _size  size of the system.
  FKPP(int _size): Size(_size)
  {
    nu=0.002;
    h=1./(Size-1.e0); uh2=nu/(h*h); rspec=4.*uh2+1.;
  }
  //! destructor.
  ~FKPP(){} 
  //! y=F(x)
  //! \param x
  //! \param y
  inline void operator()(double * Restrict x,double * Restrict  y) const
  {
    ASSUME_ALIGNED(x); ASSUME_ALIGNED(y);
    y[0]= uh2*(x[1]-x[0]) +f(x[0]);
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif 
    for(int i=1;i<Size-1;i++)
      y[i]=uh2*(x[i-1]-2.0*x[i]+x[i+1]) + f(x[i]);
    y[Size-1]= uh2*(x[Size-2]-x[Size-1]) + f(x[Size-1]);
  }
  //! return spectral radius.
  inline  double rho() const {
    return rspec;
  }
  //! initial conditions.
  //! \param x
  void init(double x[])
  {
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif     
    for(int i=0;i<Size;i++)
      x[i] = i<Size/2 ? 1.: 0.;
  }
  //! return size of the problem (number of unknowns).
  inline int size() const {return Size;}
};
#endif
