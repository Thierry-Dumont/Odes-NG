#ifndef EulerImp__h
#define EulerImp__h
#include <string>
////////////////////////////////////////////////////////////////////
/// Euler implicite.
///////////////////////////////////////////////////////////////////
namespace odes{
struct EulerImp
{
  static const int nstage=1; // number of stages.
  double A[nstage][nstage]; //A of RK method.
  double S[nstage]; //S[i]: sum of line i in A.
  double D[nstage]; // A.inv()*B, where B are the "B" of RK method.
  std::string name;
  EulerImp()
  {
    A[0][0]=1.;
    D[0]=1.;
    S[0]=1.;
    name="EulerImplicit";
  }
};
};
#endif
