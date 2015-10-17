#ifndef Rock4Coeffs__h
#define Rock4Coeffs__h
#include "coeffs-rock4"
////////////////////////////////////////////////////////////////////////
/// 
/// Wrappers to data files. Actually, this is used to make the C++
/// Rock4 code as near as possible from the fortran version.
/// We use fortran indexing.
///
////////////////////////////////////////////////////////////////////////
namespace odes {
struct Rock4Coeffs
{
  inline double recf(int i) const {return RECF[i-1];}//fortran starts at 1.
  inline double fpa(int i,int j) const 
  {
    return FPA[(i-1)*6+j-1];
  }
  inline double fpb(int i,int j) const 
  {
    return FPB[(i-1)*4+j-1];
  }
  inline double fpbe(int i,int j) const
  {
    return FPBE[(i-1)*5+j-1];
  }
  inline int ms(int i) const 
  {
    return MS[i-1];
  }
};
};
#endif
