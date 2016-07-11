#ifndef Oregonator__h
#define Oregonator__h
#include "OdesException.hpp"
#include "MatrixType.hpp"
#include "MacrosForCompilers.hpp"
///////////////////////////////////////////////////////////////////
/// Classical Oregonator problem. See Hairer & Wanner, t. 2.
//////////////////////////////////////////////////////////////////
class Oregonator
{
public:
  // there are 3 variables.
  static const int n=3;
  // the Jacobian matrix is full:
  static const int nsub=n-1;
  static const int nsup=n-1;
  // do we use Hessenberg reduction?
  static const bool Hessenberg=false;
  // compute Jacobian numerically ?
  static const bool ComputeJacobianNumerically=true;

  typedef Matrixtype<n,nsub,nsup>::Matrix Matrix;

  double u7;
  //! constructor
  Oregonator(){u7=1./77.27e0;}
  //! destructor
  ~Oregonator(){}
  //! provide a method for the initial values.
  void init(double *Restrict y)
  {
    y[0]=1.;y[1]=2.;y[2]=3.;
  }
  //! RHS for Radau5.
  inline void operator()(double t,double * Restrict y,double * Restrict res)
  {
//ASSUME_ALIGNED(y, 64);//ASSUME_ALIGNED(res, 64);
    res[0]=77.27*(y[1]+y[0]*(1.-8.375e-06*y[0]-y[1]));
    res[1]=(y[2]-(1.+y[0])*y[1])*u7;
    res[2]=0.161*(y[0]-y[2]);
  }
  inline void Jacobian(double t,const fortranVector y,const fortranVector Fy,
		       Matrix& Jac) const
  {
    Jac(1,1)=-0.00129427250*y(1) - 77.270*y(2) + 77.270;
    Jac(1,2)=-77.270*y(1) + 77.270;
    Jac(1,3)=0.0;
    Jac(2,1)=-u7*y(2);
    Jac(2,2)=-u7*(y(1)+1.0);
    Jac(2,3)=u7;
    Jac(3,1)=0.161; 
    Jac(3,2)=0.0; 
    Jac(3,3)=-0.16100;
  }
};
#endif
