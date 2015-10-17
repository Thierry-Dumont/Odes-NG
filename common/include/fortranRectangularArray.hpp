#ifndef fortranRectangularArray__h
#define fortranRectangularArray__h
#include "fortranVector.hpp"
#include "GenericException.hpp"
#include "AllocateDestroyVector.hpp"
#include <string>
#include <iostream>
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
//#define DEBUG
namespace odes{
  //////////////////////////////////////////////////////////////////////////
  /// Banded Arrray of double, of fixed size n,
  /// with kl subdiagonals and ku superdiagonals, fortran indexing, like in
  /// lapack routines (see dgbtrf for example).
  ///
  /// \brief Banded Array of double.
  //////////////////////////////////////////////////////////////////////////
  template<int n,int kl,int ku> class fortranRectangularArray
{
  static const int ldab=2*kl+ku+1;// see for example lapack/dgbtrf.
  static const int klku=kl+ku;
  static const int size=ldab*n;
  static const int l1=ldab-1;
  static const int cc=klku-ldab;
// #ifndef ON_STACK
//   double *x;
// #else
//   double x[size];
// #endif
  double *x;
public:
  //! constructor
  fortranRectangularArray()
  {
    x=allocDoubleArray(size);
// #ifndef ON_STACK
//     x=new double[size];
// #endif
  }
  //! destructor
  ~fortranRectangularArray()
  {
// #ifndef ON_STACK
//     delete[] x;
// #endif
    destroyDoubleArray(x);
  }
  //! indexing
  //! \param i
  //! \param j
  //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
  inline double operator()(int i,int j) const
  {
#ifdef DEBUG
    if(i<1||i>n||j<i-kl||j>i+ku||j<1||j>n)
      throw GenericException("fortranRectangularArray(), bad i or j:",i,j,
			     "kl=",kl,"ku=",ku);
    if(j*l1+i+cc<0||j*l1+i+cc>=size)
       throw GenericException("fortranRectangularArray(),indexing problem",
			      i,j,j*l1+i+cc);
#endif
    return *(x+j*l1+i+cc);
  }
  //! indexing.
  //! \param i
  //! \param j
  //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
  //! \note returns a reference, possible leftvalue.
  inline double& operator()(int i,int j) 
  {
#ifdef DEBUG
    if(i<1||i>n||j<i-kl||j>i+ku||j<1||j>n)
      throw GenericException("fortranRectangularArray(), bad i or j:",i,j,
			     "kl=",kl,"ku=",ku);
    if(j*l1+i+cc<0||j*l1+i+cc>=size)
       throw GenericException("fortranRectangularArray(),indexing problem",
			      i,j,j*l1+i+cc);
#endif
    return *(x+j*l1+i+cc);
  }  
  //! this =  -X.
  //! \param X the array.
  inline void equal_minus(const fortranRectangularArray<n,kl,ku>& X)
  {
#include "Ivdep.hpp"
    for(int i=0;i<size;i++)
      x[i]=-X.x[i];
  }
  //! add a value to diagonal.
  //! \param v the value.
  inline void addDiag(double v)
  {

#include "Ivdep.hpp"
    for(int i=1;i<=n;i++)
      *(x+i*ldab+cc)+=v;
  }
  //! return adress of double array.
#ifndef ON_STACK
  inline double* operator&()  {return x;}
#else
  inline double* operator&()  {return &(x[0]);}
#endif
  //! print
  void print(std::string s="")
  {
    cout<<s<<endl;
    for(int i=1;i<=n;i++)
      {
  	for(int j=MAX(1,j-kl);j<=MIN(n,j+ku);j++)
  	  cout<<operator()(i,j)<<" ";
  	cout<<endl;
      }
  }
};

}
#endif
