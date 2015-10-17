#ifndef TbbUtils__h
#define TbbUtils__h
#include "tbb/tbb.h"
#include "MacrosForCompilers.hpp"
#include "GenericException.hpp"
using namespace tbb;
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define P2(x)  ((x*x))
namespace odes{
  ////////////////////////////////////////////////////////////////////////
  /// Basic classes for Tbb: these classes are parameters to parallel_for
  /// or parallel_reduce.
  /// 1) Linear combinaisons.
  //////////////////////////////////////////////////////////////////////
  class DoLinearComb
  {
    static const unsigned int bsize=BLOCK_FACTOR/sizeof(double);
    int nvec;
    unsigned int size;
    double *v[5]; double *y;
    double coeffs[5];
    
  public:
    //! Contructor
    //! \param _size number of variables.
    //! \param cible result
    DoLinearComb(int _size,double * Restrict cible,double * Restrict  _y,
		 double a1,double * Restrict x1,
		 double a2,double * Restrict x2,
		 double a3,double * Restrict x3,
		 double a4,double * Restrict x4)
    {
      nvec=5; size=_size;
      y=_y;
      v[0]=cible; v[1]=x1; v[2]=x2; v[3]=x3,v[4]=x4;
      coeffs[1]=a1;coeffs[2]=a2; coeffs[3]=a3;coeffs[4]=a4;
    }
    DoLinearComb(int _size,double * Restrict  cible,double  * Restrict _y,
		 double a1,double  * Restrict x1,
		 double a2,double  * Restrict x2,
		 double a3,double  * Restrict x3)
    {
      nvec=4;size=_size;
      y=_y;
      v[0]=cible; v[1]=x1; v[2]=x2; v[3]=x3;
      coeffs[1]=a1;coeffs[2]=a2; coeffs[3]=a3;
    }
    DoLinearComb(int _size,double  * Restrict  cible,double * Restrict   _y,
		 double a1,double * Restrict   x1,
		 double a2,double * Restrict   x2)
    {
      nvec=3; size=_size;
      y=_y;
      v[0]=cible; v[1]=x1; v[2]=x2;
      coeffs[1]=a1;coeffs[2]=a2; 
    }
    DoLinearComb(int _size,double * Restrict cible,double* Restrict  _y,
		 double a1,double * Restrict x1)//4 parameters
    {
      nvec=2; size=_size;
      y=_y;
      v[0]=cible; v[1]=x1; 
      coeffs[1]=a1;
    }
    DoLinearComb(int _size,double  * cible, double a0,double * Restrict   _y,
		 double a1,double * Restrict   x1,
		 double a2,double *    x2)
    {
      nvec=31; size=_size;
      y=_y;
      v[0]=cible; v[1]=x1; v[2]=x2;
      coeffs[0]=a0; coeffs[1]=a1;coeffs[2]=a2; 
    }
    // copy constructor
    //! param D
    DoLinearComb(const DoLinearComb& D)
    {
      nvec=D.nvec; size=D.size;y=D.y;
      for(int i=0;i<5;i++)
	{
	  v[i]=D.v[i];
	  coeffs[i]=D.coeffs[i];
	}
    }
    //! operator used by tbb
    //! \param r a subrange of the indices.
    void operator()(const blocked_range<size_t>& r) const
    {
      int i1=r.begin()*bsize,i2=MIN(size,r.end()*bsize);
      
      double a0,a1,a2,a3,a4;

      switch(nvec)
	{
	case 1:
	  throw GenericException("DoLinearComb, nvec=",nvec);
	case 2:
	  {
	    double * Restrict x1=v[1]; a1=coeffs[1];
	    double * Restrict Y=y;
	    double * Restrict cible=v[0];
	    ASSUME_ALIGNED(x1,BLOCK_FACTOR);ASSUME_ALIGNED(cible,BLOCK_FACTOR);
	    ASSUME_ALIGNED(Y,BLOCK_FACTOR);
	    for(int i=i1;i<i2;i++)
	      cible[i]=Y[i]+ a1*x1[i];
	  }
	  break;

	case 3:
	  {
	    double * Restrict x1=v[1]; double * Restrict  x2=v[2];
	    double * Restrict Y=y; double * Restrict cible=v[0];
	    a1=coeffs[1];a2=coeffs[2];
	    ASSUME_ALIGNED(x1,BLOCK_FACTOR);ASSUME_ALIGNED(cible,BLOCK_FACTOR);
	    ASSUME_ALIGNED(Y,BLOCK_FACTOR); ASSUME_ALIGNED(x2,BLOCK_FACTOR);
	    for(int i=i1;i<i2;i++)
	      cible[i]=Y[i]+ a1*x1[i]+a2*x2[i];
	  }
	  break;	
	case 4:
	  {
	    double * Restrict x1=v[1]; double * Restrict  x2=v[2];
	    double * Restrict x3=v[3];
	    double * Restrict Y=y; double * Restrict cible=v[0];
	    ASSUME_ALIGNED(x1,BLOCK_FACTOR); ASSUME_ALIGNED(x2,BLOCK_FACTOR);
	    ASSUME_ALIGNED(x3,BLOCK_FACTOR); ASSUME_ALIGNED(cible,BLOCK_FACTOR);
	    a1=coeffs[1];a2=coeffs[2];a3=coeffs[3];
	    for(int i=i1;i<i2;i++)
	      cible[i]=Y[i]+ a1*x1[i]+a2*x2[i]+a3*x3[i];
	  }
	  break;	
	case 5:
	  {
	    double * Restrict x1=v[1]; double * Restrict x2=v[2];
	    double * Restrict x3=v[3]; double * Restrict x4=v[4];
	    double * Restrict Y=y; double * Restrict cible=v[0];
	    ASSUME_ALIGNED(x1,BLOCK_FACTOR); ASSUME_ALIGNED(x2,BLOCK_FACTOR);
	    ASSUME_ALIGNED(x3,BLOCK_FACTOR); ASSUME_ALIGNED(x4,BLOCK_FACTOR);
	    ASSUME_ALIGNED(cible,BLOCK_FACTOR);
	    a1=coeffs[1]; a2=coeffs[2]; a3=coeffs[3]; a4=coeffs[4];
	    for(int i=i1;i<i2;i++)
	      cible[i]=Y[i]+ a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
	  }
	  break;
	case 31:
	  {
	    double * Restrict x1=v[1]; double * Restrict  x2=v[2];
	    double * Restrict Y=y; double * Restrict cible=v[0];
	    a0=coeffs[0];a1=coeffs[1];a2=coeffs[2];
	    ASSUME_ALIGNED(x1,BLOCK_FACTOR);ASSUME_ALIGNED(cible,BLOCK_FACTOR);
	    ASSUME_ALIGNED(Y,BLOCK_FACTOR); ASSUME_ALIGNED(x2,BLOCK_FACTOR);
	    for(int i=i1;i<i2;i++)
	      cible[i]=a0*Y[i]+ a1*x1[i]+a2*x2[i];
	  }
	  break;	
	default:
	  throw GenericException("DoLinearComb, nvec=",nvec);
	}

    }
  };
class RPP2
{
  static const unsigned int bsize=BLOCK_FACTOR/sizeof(double);
  double *v[5]; double *y;
  double coeffs[5]; double atol,rtol; unsigned int size;
public:
  double result;
  RPP2(int _size,double _y[],double a1,double x1[],double a2,double x2[],
		double a3,double x3[],double a4,double x4[],
		double a5,double x5[],double _atol,double _rtol)
  {
    y=_y;
    v[0]=x1; v[1]=x2; v[2]=x3,v[3]=x4; v[4]=x5;
    coeffs[0]=a1;coeffs[1]=a2; coeffs[2]=a3;coeffs[3]=a4;coeffs[4]=a5;
    atol=_atol; rtol=_rtol; size=_size;
    result=0.0;
  }
 
  RPP2( RPP2& R,split):result(0)
  {
    y=R.y;
    for(int i=0;i<5;i++)
      {
	v[i]=R.v[i];
	coeffs[i]=R.coeffs[i];
      }
    atol=R.atol; rtol=R.rtol; size=R.size;
  }
  void operator()(const blocked_range<size_t>& r)
  {
    int begin=r.begin()*bsize,end=MIN(size,r.end()*bsize);
    double a1=coeffs[0],a2=coeffs[1],a3=coeffs[2],a4=coeffs[3],a5=coeffs[4];
    double * Restrict x1=v[0], * Restrict x2=v[1], * Restrict x3=v[2],
      * Restrict x4=v[3], * Restrict x5=v[4];
    ASSUME_ALIGNED(x1,BLOCK_FACTOR); ASSUME_ALIGNED(x2,BLOCK_FACTOR);
    ASSUME_ALIGNED(x3,BLOCK_FACTOR); ASSUME_ALIGNED(x4,BLOCK_FACTOR);
    ASSUME_ALIGNED(x5,BLOCK_FACTOR);

    for(int i=begin;i!=end;i++)
      {
	double ab=ABS(y[i])*rtol+atol;
	double d=1.e0/ab;
	double s=a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i]+a5*x5[i];
	result+=s*s*d*d;
      }
  }
  //! The join method for parallel_reduce
  void join(const RPP2& R)
  {
    result=result+R.result;
  }
};
};
#endif
