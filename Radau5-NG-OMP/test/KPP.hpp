#ifndef KPP__h
#define KPP__h
#include "MatrixType.hpp"
#include "MacrosForCompilers.hpp"
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
    return u*u*(1.-u);
  }
  inline double df(double u) const
  {
    return u*(2.-3.*u);
  }
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
  static const bool ComputeJacobianNumerically=true;


  typedef Matrixtype<n,nsub,nsup>::Matrix Matrix;
  //! contructor
  KPP()
  {
    //cout<<"constructor"<<endl;
    // parameters of KPP equation:
    k=1.0;
    eps=0.001;
    // 
    h=100./((double) n-1);
    eh2=eps/(h*h);
  }
  //! copy constructor
  KPP(const KPP& K)
  {
    //cout<<"copy"<<endl;
    k=K.k; eps=K.eps; h=K.h; eh2=K.eh2;
  }
  //! destructor
  ~KPP(){}
  //! provide a method for the initial condition.
  void init(double * Restrict y)
  {
    for(int i=0;i<=n/2;i++)
      y[i]=0.9999999;
    for(int i=n/2+1;i<n;i++)
      y[i]=0.0;
  }
  //! the RHS for Radau5.
  inline void operator()(double t,double * Restrict y,double * Restrict res) const
  {
    //ASSUME_ALIGNED(y, 64);//ASSUME_ALIGNED(res, 64);
    res[0]=2*eh2*(y[1]-y[0])+f(y[0]);
#include "Ivdep.hpp"
    for(int i=1;i<n-1;i++)
      res[i]=eh2*(y[i-1]-2*y[i]+y[i+1])+f(y[i]);
    res[n-1]=2*eh2*(y[n-2]-y[n-1])+f(y[n-1]);

  }
  inline void Jacobian(double t,fortranVector y,const fortranVector Fy,
		        Matrix& Jac)
  {
    Jac(1,1)=-2*eh2+df(y(1)); Jac(1,2)=2*eh2;

#include "Ivdep.hpp"
    for(int i=2;i<n;i++)
      {
	Jac(i,i)=-2.*eh2+df(y(i));
	Jac(i,i-1)=eh2; Jac(i,i+1)=eh2; 
      }
    Jac(n,n)=-2*eh2+df(y(n)); Jac(n-1,n)=2*eh2;
  }
};
#endif
