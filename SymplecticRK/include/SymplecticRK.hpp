#ifndef SymplecticRK__h
#define SymplecticRK__h
#include <cmath>
#include "GaussianMethod.hpp"
#include "Ftest.hpp"
#include "Extrap.hpp"
#include <iostream>
#define MAX(x,y) ((x)>(y)?(x):(y))
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
using namespace std;
namespace SymplectikRK{
///////////////////////////////////////////////////////////////////////////
///
/// Symplectic RK method based on Gaussian quadrature.
///
/// \brief Symplectic RK method.
//////////////////////////////////////////////////////////////////////////
template<class Fonct,int nsteps,class Double=double> class SymplecticRK
{
  static const int n=Fonct::n;
  Double v[n],Y[n*nsteps],Y1[n*nsteps],w[n],delta[n],e[n],a[n];
  int iter,itermax;
  Double Old_h;
  Fonct f;
  Extrap<nsteps,Double> Ext;
  GaussianMethod<nsteps,Double> G;
  // using GaussianMethod<nsteps,Double>::b;
  // using GaussianMethod<nsteps,Double>::A;
  // using GaussianMethod<nsteps,Double>::change_h;

  Ftest<Double>  Test;
public:
  //using GaussianMethod<nsteps,Double>::verify;
  //! constructor.
  // \param _itermax max. iteration allowed.
  SymplecticRK(int _itermax)
  { 
    itermax=_itermax;
#include "Ivdep.hpp"
    for(int i=0;i<n;i++) 
      e[i]=0.0; //for Compensated Summation.
    Old_h=1.0;
   }
  //! destructor
  ~SymplecticRK()
  {}
  //! number of iterations performed.
  inline int maxIter() const {return iter;}
  //! last computed value for the residual
  inline Double lastdiff() const {return Test.value();}
  //return a reference to F, the right hand side of the system.
  inline Fonct& rhs() {return f;}
  //! make one step.
  //! \param h time step.
  //! \param u initial value on input, value computed on output
  //! \note return true iff the fixed points iterations converged.
  bool step(Double h,Double u[])
  {
    if(Old_h!=h)
      {
	Old_h=h;
	G.change_h(h);
      }
    // this vectorizes with gcc
    for(int j=0;j<n;j++)
#include "Ivdep.hpp"
      for(int i=0;i<nsteps;i++)
	Y[i*n+j]=u[j];

    bool success=false; 

    // iteration loop (fixed point iteration):
    for(iter=0; iter<itermax;iter++)
      {
	for(int i=0;i<nsteps;i++)
	  {
#include "Ivdep.hpp"
	    for(int j=0;j<n;j++) 
	      v[j]=Y[i*n+j];
	    f(v,w);
#include "Ivdep.hpp"
	    for(int j=0;j<n;j++) Y1[i*n+j]=w[j];
	  } 
	Test=0.0;
	
	for(int i=0;i<nsteps;i++)
	  {
#include "Ivdep.hpp"
	    for(int l=0;l<n;l++)
	      v[l]=u[l];

	    // changing the order of loops vectorizes!
	    // not the classical order.
	    for(int l=0;l<n;l++)
#include "Ivdep.hpp"
	      for(int j=0;j<nsteps;j++)
		v[l]+=G.A(i,j)*Y1[j*n+l];

	    Double *Yin=Y+i*n;
#include "Ivdep.hpp"
	    for(int l=0;l<n;l++)//note: not vectorized: unhandled data-ref
	      {
		Test+=(v[l]-Yin[l])*(v[l]-Yin[l]);
	      }
#include "Ivdep.hpp"
	    for(int l=0;l<n;l++)
	      {
		Yin[l]=v[l];
	      }
	  }
	success=Test.success();
	if(success) break;
      }
   
    // 2nd phase:
    if(success)
      {
	// use Compensated Summation (see Hairer & Co. 2nd edition, page 323).
	for(int i=0;i<n;i++)
	  delta[i]=0.0;
	
	  for(int i=0;i<nsteps;i++)
	    {
	      double bi=G.b[i];
#include "Ivdep.hpp"
	      for(int j=0;j<n;j++)
		delta[j]+=bi*Y1[i*n+j];
	    }
#include "Ivdep.hpp"
	for(int i=0;i<n;i++)
	  a[i]=u[i];
#include "Ivdep.hpp"
	for(int i=0;i<n;i++)
	  e[i]+=delta[i];
#include "Ivdep.hpp"
	for(int i=0;i<n;i++)
	  u[i]=a[i]+e[i];
#include "Ivdep.hpp"
	for(int i=0;i<n;i++)
	  e[i]+=(a[i]-u[i]);
      }
    return success;
  }


};
}
#endif
