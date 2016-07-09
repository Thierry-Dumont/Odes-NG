#ifndef RKMethod2__h
#define RKMethod2__h
#include <string>
//////////////////////////////////////////////////////////////////////
///  Crouzeix method. cf. Hairer-Wanner page 100, figure 6.18
//////////////////////////////////////////////////////////////////////
namespace odes{
struct RKMethod2
{
  static const int nstage=3; //!< number of stages.
  double A[nstage][nstage]; //!< A of RK method.
  double S[nstage]; //!< S[i]: sum of line i in A.
  double D[nstage]; //!< A.inv()*B, where B are the "B" of RK method.
  std::string name;
  RKMethod2()
  {
    A[0][0]=1.06857902130163;
    A[1][0]=-0.568579021301629;
    A[1][1]=1.06857902130163;
    A[2][0]=2.13715804260326;
    A[2][1]=-3.27431608520651;
    A[2][2]=1.06857902130163;

    D[0]=0.445622407287714;
    D[1]=1.06417777247591;
    D[2]=0.120614758428183;

    S[0]=1.06857902130163;
    S[1]=0.500000000000000;
    S[2]=-0.0685790213016287;

    name="CrouzeixRaviart-A";

  }
};
};
#endif
