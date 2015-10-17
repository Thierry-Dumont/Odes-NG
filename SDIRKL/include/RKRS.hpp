#ifndef RKRS__h
#define RKRS__h
#include <string>
/////////////////////////////////////////////////////////////////////////////
/// Méthode of Ropp and Shadid JCP 203 (2005) 449-466.
/////////////////////////////////////////////////////////////////////////////
struct RKRS
{
  static const int nstage=2; //!< number of stages.
  double A[nstage][nstage]; //!< \f$ A_{i,j}\f$ coefficients of the RK method.
  double S[nstage]; //!< S[i]: sum of line i in A: \f$ S_i= \sum_j A_{i,j}\f$
  double D[nstage]; //!< \f$ D= B.A^{-1} \f$ where B are the "B" of RK method.
  std::string name;
  RKRS()
  {
    A[0][0]=1.70710678118655;
    A[1][0]=-2.41421356237309;
    A[1][1]=1.70710678118655;

    D[0]=0.707106781186547;
    D[1]=0.292893218813452;
 
    S[0]=1.70710678118655;
    S[1]=-0.707106781186548;

    name="RoppShadid";
    
  }
};
#endif
