#ifndef MatrixType__h
#define MatrixType__h
#include "Matrices.hpp"
/////////////////////////////////////////////////////////////////////////////
/// Helper class, used at compile time to determine the type of Jacobian
/// matrices.
/// \brief   determine the type of Jacobian matrices.
////////////////////////////////////////////////////////////////////////////
namespace odes
{
template<int n,int nsub,int nsup> struct Matrixtype
{
  static const bool Full=(n-nsub)==1&&(n-nsup)==1;
  typedef typename Matrices<Full,false,n,nsub,nsup>::MatrixReal Matrix;
};
}
#endif
