#ifndef KPP__h
#define KPP__h
#include "MatrixType.hpp"
/////////////////////////////////////////////////////////////////////////
/// Kolmogorov Piskunov Petrovskii réaction diffusion in 1d.
/// u_t = \eps u_xx +  k u^2 (1-u).
/// with Neuman Boundary conditions. Finite differences discretization.
///////////////////////////////////////////////////////////////////////
class KPP
{
  double k,eps,h,eh2;
  inline double f(double u) const
  {
    return k*u*u*(1.-u);
  }
  inline double df(double u) const
  {
    return k*u*(2.-3.*u);
  }
public:
  // n: number of points in the finite difference discretization 
  //    this is also the number of variables.
  static const int n=100000; 
  // the Jacobian is tridiagonal:
  static const int nsub=1;
  static const int nsup=1;
  // do we use Hessenberg reduction? (we cannot, because Jacobian matrix is not
  // full).
  static const bool Hessenberg=false;
  // compute Jacobian numerically ?
  static const bool ComputeJacobianNumerically=false;
  //is the system autonomous ?
  static const bool autonomous=true;
  //if the systeme is not autonomous, do we provide the derivative of F
  //with respect to t ? Otherwise it will be computed numerically.
  static const bool use_DF_t=false;
  // which method do we use in Rodas?
  static const int method=1;

  typedef Matrixtype<n,nsub,nsup>::Matrix Matrix;
  //! contructor
  KPP()
  {
    // parameters of KPP equation:
    k=1.0;
    eps=0.001;
    // 
    h=100./((double) n-1);
    eh2=eps/(h*h);
  }
  //! destructor
  ~KPP(){}
  //! provide a method for the initial condition.
  void init(double y[])
  {
    for(int i=0;i<=n/2;i++)
      y[i]=0.9999999;
    for(int i=n/2+1;i<n;i++)
      y[i]=0.0;
  }
  //! the RHS for Radau5.
  inline void operator()(double t,double y[],double res[]) const
  {
    res[0]=2*eh2*(y[1]-y[0])+f(y[0]);
    for(int i=1;i<n-1;i++)
      res[i]=eh2*(y[i-1]-2*y[i]+y[i+1])+f(y[i]);
    res[n-1]=2*eh2*(y[n-2]-y[n-1])+f(y[n-1]);

  }
  inline void Jacobian(double t,fortranVector y,const fortranVector Fy,
		        Matrix& Jac)
  {
    Jac(1,1)=-2*eh2+df(y(1)); Jac(1,2)=2*eh2;
    for(int i=2;i<n;i++)
      {
	Jac(i,i)=-2.*eh2+df(y(i));
	Jac(i,i-1)=eh2; Jac(i,i+1)=eh2; 
      }
    Jac(n,n)=-2*eh2+df(y(n)); Jac(n,n)=2*eh2;
  }
  inline void DF_t(double t,const fortranVector y,fortranVector& DFt)
  {
    throw GenericException("Oregonator:  DF_t called");
  }
};
#endif
