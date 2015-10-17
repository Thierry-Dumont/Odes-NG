#ifndef FRadau5__h
#define FRadau5__h
#include "MatrixType.hpp"
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include "binit.hpp"
/////////////////////////////////////////////////////////////////////////
/// u_t = \eps u_xx +  f
/// with Neuman Boundary conditions. Finite differences discretization.
///////////////////////////////////////////////////////////////////////
class FRadau5
{
  double k,eps,h,eh2;
  double *b;
public:
  // n: number of points in the finite difference discretization 
  //    this is also the number of variables.
  static const int n=100; 
  // the Jacobian is tridiagonal:
  static const int nsub=1;
  static const int nsup=1;
  // do we use Hessenberg reduction? (we cannot, because Jacobian matrix is not
  // full).
  static const bool Hessenberg=false;
  // compute Jacobian numerically ?
  static const bool ComputeJacobianNumerically=false;


  typedef Matrixtype<n,nsub,nsup>::Matrix Matrix;
  //! contructor
  FRadau5()
  {
    // parameters of FRadau5 equation:
    k=1.0;
    eps=1.;
    // 
    h=1./((double) n-1);
    eh2=eps/(h*h);
    b=allocDoubleArray(n);
    binit(b,n);
  }
  //! destructor
  ~FRadau5()
  {
    destroyDoubleArray(b);
  }
  //! provide a method for the initial condition.
  void init(double * Restrict y)
  {
    cout<<"dans init"<<endl;
    for(int i=0;i<=n/2;i++)
      y[i]=0.9999999;
    for(int i=n/2+1;i<n;i++)
      y[i]=0.0;
  }
  //! the RHS for Radau5.
  inline void operator()(double t,double * Restrict y,double * Restrict res) const
  {
    ASSUME_ALIGNED(y);ASSUME_ALIGNED(res);
    res[0]=2*eh2*(y[1]-y[0])+b[0];
#include "Ivdep.hpp"
 for(int i=1;i<n-1;i++)
   res[i]=eh2*(y[i-1]-2*y[i]+y[i+1])+b[i];
 res[n-1]=2*eh2*(y[n-2]-y[n-1])+b[n-1];

  }
  inline void Jacobian(double t,fortranVector y,const fortranVector Fy,
		        Matrix& Jac)
  {
    Jac(1,1)=-2*eh2; Jac(1,2)=2*eh2;
#include "Ivdep.hpp"
    for(int i=2;i<n;i++)
      {
	Jac(i,i)=-2.*eh2;
	Jac(i,i-1)=eh2; Jac(i,i+1)=eh2; 
      }
    Jac(n,n)=-2*eh2; Jac(n-1,n)=2*eh2;
  }
};
#endif
