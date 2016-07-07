#ifndef fortranRectangularComplexArray__h
#define fortranRectangularComplexArray__h
#include "fortranVector.hpp"
#include "OdesException.hpp"
#include "AllocateDestroyVector.hpp"
#include <string>
#include <iostream>
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
//#define DEBUG
namespace odes{
  ////////////////////////////////////////////////////////////////////////
  /// Banded Array of Complex, of fixed size n,
  /// with kl subdiagonals and ku superdiagonals, fortran indexing, like in
  /// lapack routines (see zgbtrf for axample).
  /// Actually, we do not make use of complex<double> class (we just want
  /// to store coefficients for lapack zge* routines).
  /// \brief Banded Array of complex..
  ////////////////////////////////////////////////////////////////////////
  template<int n,int kl,int ku> class fortranRectangularComplexArray
  {
    static const int ldab=2*kl+ku+1;// see for example lapack/dgbtrf.
    static const int klku=kl+ku;
    static const int size=2*ldab*n;
    static const int l1=ldab-1;
    static const int cc=klku-ldab;
// #ifndef ON_STACK
//   double *x;
// #else
//   double x[size];
// #endif
    double *x;
  public:
    //! contructor
    fortranRectangularComplexArray()
    {
// #ifndef ON_STACK
//     x=new double[size];
// #endif
      x=allocDoubleArray(size);
    }
    //! destructor
    ~fortranRectangularComplexArray()
    {
// #ifndef ON_STACK
//     delete[] x;
// #endif
      destroyDoubleArray(x);
    }
    //set Real and Imaginary part, for (i,j):
    inline void set(int i,int j,double RealPart,double ImagPart=0)
    {
#ifdef DEBUG
      if(i<1||i>n||j<i-kl||j>i+ku||j<1||j>n)
	throw OdesException("fortranRectangularComplexArray(), bad i or j:"
			       ,i,j,"kl=",kl,"ku=",ku);
      if(j*l1+i+cc<0||2*(j*l1+i+cc)>=size-1)
	throw OdesException("fortranRectangularComplexArray()",
			       "indexing problem",i,j,2*(j*l1+i+cc),size);
#endif
    
      int k=2*(j*l1+i+cc);
      x[k]=RealPart; x[k+1]=ImagPart;
    }
    //!Indexing, return a reference to "real" part.
    //! \param i
    //! \param j
    //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
    inline double& Re(int i,int j)
    {
#ifdef DEBUG
    if(i<1||i>n||j<i-kl||j>i+ku||j<1||j>n)
      throw OdesException("fortranRectangularComplexArray(), bad i or j:"
			     ,i,j,"kl=",kl,"ku=",ku);
    if(j*l1+i+cc<0||2*(j*l1+i+cc)>=size-1)
       throw OdesException("fortranRectangularComplexArray()",
			      "indexing problem",i,j,j*l1+i+cc,size);
#endif
       return x[2*(j*l1+i+cc)];
    }
    //!Indexing, return a reference to "Imaginary" part.
    //! \param i
    //! \param j
    //! \note fortran indexing (row major, i>=1, j>=1,i<=n, j<=n).
    inline double& Im(int i,int j)
    {
#ifdef DEBUG
    if(i<1||i>n||j<i-kl||j>i+ku||j<1||j>n)
      throw OdesException("fortranRectangularComplexArray(), bad i or j:"
			     ,i,j,"kl=",kl,"ku=",ku);
    if(j*l1+i+cc<0||2*(j*l1+i+cc)>=size-1)
       throw OdesException("fortranRectangularComplexArray()",
			      "indexing problem",i,j,j*l1+i+cc);
#endif
       return x[2*(j*l1+i+cc)+1];
    }
    //! return adress of double array.
#ifndef ON_STACK
  inline double* operator&()  {return x;}
#else
  inline double* operator&()  {return &(x[0]);}
#endif
    //!print
    void print(std::string s="")
    {
      cout<<s<<endl;
      for(int i=1;i<=n;i++)
    	{
    	  for(int j=MAX(1,j-kl);j<=MIN(n,j+ku);j++)
    	    cout<<Re(i,j)<<" ";
    	  cout<<endl;
    	}
      cout<<"---"<<endl;
      for(int i=1;i<=n;i++)
    	{
    	  for(int j=MAX(1,j-kl);j<=MIN(n,j+ku);j++)
    	    cout<<Im(i,j)<<" ";
    	  cout<<endl;
    	}    
    } 
  };
}
#endif
