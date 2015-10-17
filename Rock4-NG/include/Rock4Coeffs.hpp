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
  inline double fpa(int i,int j) const //fortran indexing
  {
    return FPA[(j-1)*50+i-1];
  }
  inline double fpb(int i,int j) const //fortran indexing
  {
    return FPB[(j-1)*50+i-1];
  }
  inline double fpbe(int i,int j) const //fortran indexing
  {
    return FPBE[(j-1)*50+i-1];
  }
  inline int ms(int i) const //fortran indexing
  {
    return MS[i-1];
  }
};
};
#endif
