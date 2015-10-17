#ifndef fortranArray__h
#define fortranArray__h
#include "fortranVector.hpp"
#include "GenericException.hpp"
#include <string>
#include <iostream>
#include "AllocateDestroyVector.hpp"
namespace odes{
  ////////////////////////////////////////////////////////////////////////
  /// Array of double, of fixed size n.n, fortran indexing.
  ///
  /// \brief Array of double.
  ///
  ////////////////////////////////////////////////////////////////////////

  template<int n> class fortranArray
{
  static const int n2=n*n;
  double *x;
public:
  //! constructor
  fortranArray()
  {
    x=allocDoubleArray(n2);
  }
  //! constructor (copy).
  fortranArray(const fortranArray<n>& AA )
  {
    ASSUME_ALIGNED(x);ASSUME_ALIGNED(AA.x);
#include "Ivdep.hpp"
    for(int i=0;i<n*n;i++)
      x[i]=AA.x[i];
  }
  //! destructor
  ~fortranArray()
  {
    destroyDoubleArray(x);
  }
  //! indexing
  //! \param i
  //! \param j
  //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
  inline double operator()(int i,int j) const
  {
#ifdef DEBUG
    if(i<1||i>n||j<1||j>n)
      throw GenericException("fortranArray(), bad i or j",i,j);
#endif
    return *(x+(j-1)*n+i-1);
  }
  //! indexing.
  //! \param i
  //! \param j
  //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
  //! \note returns a reference, possible leftvalue.
  inline double& operator()(int i,int j) 
  {
#ifdef DEBUG
    if(i<1||i>n||j<1||j>n)
      throw GenericException("fortranArray(&), bad i or j",i,j);
#endif
    return *(x+(j-1)*n+i-1);
  }  
  //! this = -array.
  //! \param X the array.
  inline void equal_minus(fortranArray<n>& X)
  {
    ASSUME_ALIGNED(x);ASSUME_ALIGNED(X.x);
#include "Ivdep.hpp"
    for(int i=0;i<n*n;i++)
      x[i]=-X.x[i];
  }
  //! add a value to diagonal.
  //! \param v the value.
  inline void addDiag(double v)
  {
#include "Ivdep.hpp"
    for(int i=0;i<n;i++)
      *(x+i*(n+1))+=v;
  }
  //! return adress of double array.
#ifndef ON_STACK
  inline double* operator&()  {return x;}
#else
  inline double* operator&()  {return &(x[0]);}
#endif

  //! print
  void print(std::string s="") const
  {
    cout<<s<<endl;
    for(int i=1;i<=n;i++)
      {
	for(int j=1;j<=n;j++)
	  cout<<operator()(i,j)<<" ";
	cout<<endl;
      }
  }
};

}
#endif
