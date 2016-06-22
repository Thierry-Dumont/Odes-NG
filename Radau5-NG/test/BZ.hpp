#ifndef BZ__h
#define BZ__h
#include "MatrixType.hpp"
#include "MacrosForCompilers.hpp"
/////////////////////////////////////////////////////////////////////////
/// BZ system of 3 reaction diffusion equations in dmension 1.
/// with Neuman Boundary conditions. Finite differences discretization.
///////////////////////////////////////////////////////////////////////
class BZ
{
  typedef Matrixtype<3,2,2>::Matrix Matrix3;
  double k,h;
  double u7;
  double eps[3];
  inline double f1(double u,double v,double w) const
  {
    return 77.27*(v+u*(1.-8.375e-06*u-v));
  }
  inline double f2(double u,double v,double w) const
  {
    return (w-(1.+u)*v)*u7;
  }
  inline double f3(double u,double v,double w) const
  {
    return 0.161*(u-w);
  }
  inline void df(const fortranVector& y,Matrix3& Jac) const
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
public:
  // nd: number of points in the finite difference discretization 
  static const int nd=100;
  // n:  the number of variables.
  static const int n=3*nd; 
  // the Jacobian is full:
  static const int nsub=2;
  static const int nsup=2;
  // do we use Hessenberg reduction? 
  static const bool Hessenberg=true;
  // compute Jacobian numerically ?
  static const bool ComputeJacobianNumerically=false;


  typedef Matrixtype<n,nsub,nsup>::Matrix Matrix;
 

  //! contructor
  BZ()
  {
    // 
    h=1./((double) nd-1);
    u7=1./77.27e0;
    //diffusion
    eps[0]=0.0025/(h*h); eps[1]=eps[0]; eps[2]=0.0015/(h*h);
  }
  //! destructor
  ~BZ(){}
  //! provide a method for the initial condition.
  void init(double * Restrict y)
  {
    for(int i=0;i<n;i++)
      y[i]=0.0;
    for(int i=nd/2-nd/6;i<nd/2+nd/6;i++)
      y[i]=1.0;
  }
  //! the RHS for Radau5.
  inline void operator()(double t,double * Restrict y,
			 double * Restrict res) const
  {
    double * Restrict res1=res+nd;
    double * Restrict y1=y+nd;
    double * Restrict res2=res1+nd;
    double * Restrict y2=y1+nd;
    //
    double eh2=eps[0];
    ASSUME_ALIGNED(y);ASSUME_ALIGNED(res);
    res[0]=2*eh2*(y[1]-y[0])+f1(y[0],y1[0],y2[0]);
#include "Ivdep.hpp"
    for(int i=1;i<nd-1;i++)
      res[i]=eh2*(y[i-1]-2*y[i]+y[i+1])+f1(y[i],y1[i],y2[i]);
    res[nd-1]=2*eh2*(y[nd-2]-y[nd-1])+f1(y[nd-1],y1[nd-1],y2[nd-1]);
    //
    eh2=eps[1];
    res1[0]=2*eh2*(y1[1]-y1[0])+f2(y[0],y1[0],y2[0]);
#include "Ivdep.hpp"
    for(int i=1;i<nd-1;i++)
      res1[i]=eh2*(y1[i-1]-2*y1[i]+y1[i+1])+f2(y[i],y1[i],y2[i]);
    res1[nd-1]=2*eh2*(y1[nd-2]-y1[nd-1])+f2(y[nd-1],y1[nd-1],y2[nd-1]);
    //
    eh2=eps[2];
    res2[0]=2*eh2*(y2[1]-y2[0])+f3(y[0],y1[0],y2[0]);
#include "Ivdep.hpp"
    for(int i=1;i<nd-1;i++)
      res2[i]=eh2*(y2[i-1]-2*y2[i]+y2[i+1])+f3(y[i],y1[i],y2[i]);
    res2[nd-1]=2*eh2*(y2[nd-2]-y2[nd-1])+f3(y[nd-1],y1[nd-1],y2[nd-1]);
  }
  inline void Jacobian(double t,fortranVector y,const fortranVector Fy,
		        Matrix& Jac)
  {
    // for(int j=1;j<=n;j++)
    //   for(int i=1;i<=n;i++)
    // 	Jac(i,j)=0.0;
    //linear part.
    double eh2=eps[0];
    Jac(1,1)=-2*eh2; Jac(1,2)=2*eh2;
#include "Ivdep.hpp"
    for(int i=2;i<nd;i++)
      {
	Jac(i,i)=-2.*eh2;
	Jac(i,i-1)=eh2; Jac(i,i+1)=eh2; 
      }
    Jac(nd,nd)=-2*eh2; Jac(nd-1,nd)=2*eh2;

    //
    eh2=eps[1];
    Jac(nd+1,nd+1)=-2*eh2; Jac(nd+1,nd+2)=2*eh2;
    for(int i=2;i<nd;i++)
      {
	Jac(nd+i,nd+i)=-2.*eh2;
	Jac(nd+i,nd+i-1)=eh2; Jac(nd+i,nd+i+1)=eh2; 
      }
    Jac(2*nd,2*nd)=-2*eh2; Jac(2*nd-1,2*nd)=2*eh2;
    //
    eh2=eps[1];
    Jac(2*nd+1,2*nd+1)=-2*eh2; Jac(2*nd+1,2*nd+2)=2*eh2;
    for(int i=2;i<nd;i++)
      {
	Jac(2*nd+i,2*nd+i)=-2.*eh2;
	Jac(2*nd+i,2*nd+i-1)=eh2; Jac(2*nd+i,2*nd+i+1)=eh2; 
      }
    Jac(3*nd,3*nd)=-2*eh2; Jac(3*nd-1,3*nd)=2*eh2;

    // non linear part:
    Matrix3 Jac3;fortranVector yv(3);
    
    for(int i=1;i<=nd;i++)
      {
	yv(1)=y(i);yv(2)=y(i+nd);yv(3)=y(i+2*nd);
	df(yv,Jac3);
	for(int j1=1;j1<=3;j1++)
	  for(int i1=1;i1<=3;i1++)
	    Jac(i+(i1-1)*nd,i+(j1-1)*nd)+=Jac3(i1,j1);
      }
  }
};
#endif
