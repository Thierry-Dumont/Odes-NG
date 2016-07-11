#ifndef E5__h
#define E5__h
#include "MatrixType.hpp"
///////////////////////////////////////////////////////////////////
/// Classical E5 problem. See Hairer & Wanner, t. 2.
//////////////////////////////////////////////////////////////////
class E5
{
public:
  // there are 3 variables.
  static const int n=4;
  // the Jacobian matrix is full:
  static const int nsub=n-1;
  static const int nsup=n-1;
  // do we use Hessenberg reduction?
  static const bool Hessenberg=true;
  // compute Jacobian numerically ?
  static const bool ComputeJacobianNumerically=true;

  typedef Matrixtype<n,nsub,nsup>::Matrix Matrix;
  // constants:
  const double M=1.e+06;
  const double C=1.13e+03;
  const double A=7.89e-10;
  const double B=1.1e+07;
  double MC;// g++ does not like MC as static const double MC=M*C.
  //! constructor
  E5()
  {
    MC=M*C;
  }
  //! destructor
  ~E5(){}
  //! provide a method for the initial values.
  void init(double y[])
  {
    y[0]=1.76e-03,y[1]=0.0;y[2]=0.0;y[3]=0.0;
  }
  //! RHS for Radau5.
  inline void operator()(double t,double y[],double res[]) const
  {
    res[0]=-A*y[0]-B*y[0]*y[2];
    res[1]= A*y[0]            -MC*y[1]*y[2];
    res[3]=        B*y[0]*y[2]             -C*y[3];
    // do not compute res[2] like this:
    //res[2]= A*y[0]-B*y[0]*y[2]-MC*y[1]*y[2]+C*y[3];
    //but use use the invariant (see Hairer & Wanner t2):

    res[2]=res[1]-res[3];
  }
  inline void Jacobian(double t,fortranVector y,const fortranVector Fy,
		       Matrix& Jac)
  {
    // fake Jacobian.
    throw OdesException("Analytical Jacobian not coded");
  }
};
#endif
