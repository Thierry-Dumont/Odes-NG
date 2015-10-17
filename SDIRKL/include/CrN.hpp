#ifndef CrN__h
#define CrN__h
#include <string>
/////////////////////////////////////////////////////////////////////////////
/// Crank-Nicolson method.
/////////////////////////////////////////////////////////////////////////////
struct CrN
{
  static const int nstage=1; // number of stages.
  double A[nstage][nstage]; //A of RK method.
  double S[nstage]; //S[i]: sum of line i in A.
  double D[nstage]; // A.inv()*B, where B are the "B" of RK method.
  std::string name;
  CrN()
  {
    A[0][0]=0.5;
    D[0]=1.;
    S[0]=0.5;
    name="CrN";
  }
};
#endif
