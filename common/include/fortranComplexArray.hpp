#ifndef fortranComplexArray__h
#define fortranComplexArray__h
#include "fortranVector.hpp"
#include "OdesException.hpp"
#include <string>
#include <iostream>
#include "AllocateDestroyVector.hpp"
namespace odes{
  ////////////////////////////////////////////////////////////////////////
  /// Array of Complex, of fixed size n.n, fortran indexing.
  /// Actually, we do not make use of complex<double> class (we just want
  /// to store coefficients for lapack zge* routines).
  /// \brief Array of complex..
  ////////////////////////////////////////////////////////////////////////
  template<int n> class fortranComplexArray
  {
// #ifndef ON_STACK
//     double *x;
// #else
//     double x[2*n*n];
// #endif
    double *x;
    static const int d2=2*n;
    static const int size=d2*n;
  public:
    //! contructor
    fortranComplexArray()
    {
// #ifndef ON_STACK
//     x=new double[size];
// #endif
      x=allocDoubleArray(size);
    }
    //! destructor
    ~fortranComplexArray()
    {
// #ifndef ON_STACK
//     delete[] x;
// #endif
      destroyDoubleArray(x);
    }
    //set Real and Imaginary part, for (i,j):
    inline void set(int i,int j,double RealPart,double ImagPart=0)
    {
      int k=(j-1)*d2+2*(i-1);
      x[k]=RealPart; x[k+1]=ImagPart;
    }
    //!Indexing, return a reference to "real" part.
    //! \param i
    //! \param j
    //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
    inline double& Re(int i,int j)
    {
      return x[(j-1)*d2+2*(i-1)];
    }
    //!Indexing, return a reference to "Imaginary" part.
    //! \param i
    //! \param j
    //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
    inline double& Im(int i,int j)
    {
      return x[(j-1)*d2+2*i-1];
    }
    //! return adress of double array.
 
#ifndef ON_STACK
  inline double* operator&()  {return x;}
#else
  inline double* operator&()  {return &(x[0]);}
#endif  //!print
    void print(std::string s="")
    {
      cout<<s<<endl;
      for(int i=1;i<=n;i++)
	{
	  for(int j=1;j<=n;j++)
	    cout<<Re(i,j)<<" ";
	  cout<<endl;
	}
      cout<<"---"<<endl;
      for(int i=1;i<=n;i++)
	{
	  for(int j=1;j<=n;j++)
	    cout<<Im(i,j)<<" ";
	  cout<<endl;
	}    
    } 
  };

}
#endif
