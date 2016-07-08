#ifndef GaussianMethod__h
#define GaussianMethod__h
#include <cmath>
#include <iostream>
#include "GenericException.hpp"
#include "icoeffs.hpp"
//#define DEBUG
namespace SymplectikRK{
  ///////////////////////////////////////////////////////////////////////////
  /// Build coefficients for symplectic gaussian RK.
  /// \brief  Build coefficients for symplectic gaussian RK.
  ///////////////////////////////////////////////////////////////////////////
  using namespace odes;
  template<int nsteps,class Double> class GaussianMethod
  {
    // default template, never instanciated.
  };
  // Classical "double" precision. 
  template<int nsteps> class GaussianMethod<nsteps,double>
  {
  protected:
    double a[nsteps*nsteps];
    double b[nsteps];//,d[nsteps];
  
    //! constructor.
    GaussianMethod()
    {
      //init();
      icoeffs<nsteps,double>(a,b);
    }
    //! multiply coefficients by h.
    //! \param h
    void change_h(double h)
    {
      icoeffs<nsteps,double>(a,b);
      for(int i=0;i<nsteps*nsteps;i++)
	a[i]*=h;
      for(int i=0;i<nsteps;i++)
	b[i]*=h;
    }
    
    //! destructor.
    ~GaussianMethod(){}
    //! access to an array.
    inline double A(int i,int j) const
    {
#ifdef DEBUG
      if(i<0||i>=nsteps||j<0||j>nsteps)
	throw GenericException("SymplecticRK::GaussianMethod::A(int i,int i), i",
			       i," > nsteps or < 0");
#endif
      return a[i*nsteps+j];
    }
  public:
    //! verify theorem of Cooper (see Hairer & Co)
    void verify()
    {
      for(int i=0;i<nsteps;i++)
	{
	  for(int j=0;j<nsteps;j++)
	    std::cout<<b[i]*A(i,j)+b[j]*A(j,i)-b[i]*b[j]<<" ";
	  std::cout<<std::endl;
	}
    }
  };
  // long double of g++ on X86/64.
  template<int nsteps> class GaussianMethod<nsteps,long double>
  {
  protected:
    typedef long double Double;
    Double a[nsteps*nsteps];
    Double b[nsteps],d[nsteps];
  
    //! constructor.
    GaussianMethod()
    {
      //init();
      icoeffs<nsteps,long double>(a,b);
    }
    //! multiply coefficients by h.
    //! \param h
    void change_h(double h)
    {
      //init();
      icoeffs<nsteps,long double>(a,b);
      for(int i=0;i<nsteps*nsteps;i++)
	a[i]*=h;
      for(int i=0;i<nsteps;i++)
	b[i]*=h;
    }

    //! destructor.
    ~GaussianMethod(){}
    //! access to a array.
    inline Double A(int i,int j) const
    {
#ifdef DEBUG
      if(i<0||i>=nsteps||j<0||j>nsteps)
	throw GenericException("SymplecticRK::GaussianMethod::A(int i,int i), i",
			       i," > nsteps or < 0");
#endif
      return a[i*nsteps+j];
    }
  public:
    //! verify theorem of Cooper (see Hairer & Co)
    void verify()
    {
      for(int i=0;i<nsteps;i++)
	{
	  for(int j=0;j<nsteps;j++)
	    std::cout<<b[i]*A(i,j)+b[j]*A(j,i)-b[i]*b[j]<<" ";
	  std::cout<<std::endl;
	}
    }
  };
}
#endif
