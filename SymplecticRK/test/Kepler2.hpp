#ifndef Kepler2__h
#define Kepler2__h
#include <cmath>
template<class Double> class Kepler2
{
 public:
  static const int n=4;
  Kepler2()
  {}
  inline void operator()(Double X[],Double Y[]) const
  {
    //symplify notations with references:
    Double& q1=X[0]; Double& q2=X[1];
    Double& p1=X[2]; Double& p2=X[3];

    Y[2]=+q1/sqrt(pow(q1*q1 + q2*q2,3));
    Y[3]=+q2/sqrt(pow(q1*q1 + q2*q2,3));
    Y[0]=-p1;
    Y[1]=-p2;
  }
  inline Double H(Double X[]) const //Hairer and Co. page 9..
  {
    //symplify notations with references:
    Double& q1=X[0]; Double& q2=X[1];
    Double& p1=X[2]; Double& p2=X[3];

    return 0.5*(p1*p1+p2*p2)-1./sqrt(q1*q1+q2*q2);//Hamiltonian.
    //return q1*p2-q2*p1;// momentum.
  }
  inline void init(Double U[])
  {
    U[0]=0.1;U[1]=0.; U[2]=0.05;U[3]=0.1;
  }
};
#endif
