#pragma once
#include <iostream>
#include "Rock4Coeffs.hpp"
#include "MacrosForCompilers.hpp"
#include "OdesException.hpp"
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
using namespace std;
namespace odes {
  ///////////////////////////////////////////////////////////////////////////
  /// Perform basic operations for Rock4. This is the both for sequential
  /// and openmp computations.
  ///
  //////////////////////////////////////////////////////////////////////////
  template<class Fonc> class ToolOMP
  {
    int mr,mz,size; double h;
    Rock4Coeffs Coeffs;
  public:
    void print(double x[],string s="")
    {
      std::cout.precision(16);
      cout<<endl<<s<<" ";
      for(int i=0;i<size;i++)
	cout<<scientific<<x[i]<<" ";
      cout<<endl;
    }
    //! Constructor
    //! \param _mr see Rock4.hpp
    //! \param _mz see Rock4.hpp
    //! \param _size see Rock4.hpp
    //! \param h time step.
    ToolOMP(int _mr,int _mz,int _size,double _h):
      mr(_mr),mz(_mz),size(_size),h(_h){}
    //! Destructor
    ~ToolOMP(){}
    //! The recurrence for one degree.
    //! \param F RHS.
    //! \param g1 preceeding term
    //! \param g0 on input: 2 terms ago; result on output.
    //! \param deg degree..
    void apply(Fonc& F,int deg,double * Restrict g1, double * Restrict g0,
	       double * Restrict vtemp)
    {
      double temp1=h*Coeffs.recf(mr+2*(deg-2)+1),
	temp3=-Coeffs.recf(mr+2*(deg-2)+2);
      double temp2=1.0-temp3;
      //
      F(g1,vtemp);
      ASSUME_ALIGNED(g0); ASSUME_ALIGNED(g1);
      ASSUME_ALIGNED(vtemp);
      //VEC!
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif      
      for(int i=0;i<size;i++)
	g0[i]=temp1*vtemp[i]+temp2*g1[i]+temp3*g0[i];
    }
    void apply(Fonc& F,double * Restrict g2, double * Restrict g1,
	       double * Restrict g0)
    {
      // here, degree==2.
      double temp1=h*Coeffs.recf(mr+1),	temp3=-Coeffs.recf(mr+2);
      double temp2=1.0-temp3;

      F(g1,g2);
      ASSUME_ALIGNED(g0); ASSUME_ALIGNED(g1);
      ASSUME_ALIGNED(g2);
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif 
      for(int i=0;i<size;i++)
	g2[i]=temp1*g2[i]+temp2*g1[i]+temp3*g0[i];
    }
    void cl(double * Restrict g,double * Restrict y,double * Restrict fn)
    {
      double temp=h*Coeffs.recf(mr);
      ASSUME_ALIGNED(g);  ASSUME_ALIGNED(y);
      ASSUME_ALIGNED(fn);
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif 
      for(int i=0;i<size;i++)
	g[i]=y[i]+temp*fn[i];
    }
    void cl(double * Restrict cible,double * Restrict y,double a1,
	    double *Restrict x1)
    {
      ASSUME_ALIGNED(cible);  ASSUME_ALIGNED(y);
      ASSUME_ALIGNED(x1);
      //VEC!
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i];
    }
    void cl(double  * Restrict cible,double  * Restrict y,
	    double a1,double  * Restrict x1,
	    double a2,double  * Restrict x2)
    {
      ASSUME_ALIGNED(cible);  ASSUME_ALIGNED(y);  
      ASSUME_ALIGNED(x1);  ASSUME_ALIGNED(x2); 
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif 
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i]+a2*x2[i];
    }
    void cl(double * Restrict cible,double * Restrict y,
	    double a1,double* Restrict  x1,
	    double a2,double * Restrict x2,double a3,double * Restrict x3)
    {
      ASSUME_ALIGNED(cible);  ASSUME_ALIGNED(y);  
      ASSUME_ALIGNED(x1);  ASSUME_ALIGNED(x2);
      ASSUME_ALIGNED(x3); 
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif 
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i];
    }
    void ct(double * Restrict  cible,double * Restrict  y,
	    double a1,double * Restrict  x1,double a2,double * Restrict  x2,
	    double a3,double * Restrict  x3,double a4,double  * Restrict x4)
    {
      ASSUME_ALIGNED(cible);  ASSUME_ALIGNED(y);  
      ASSUME_ALIGNED(x1);  ASSUME_ALIGNED(x2);
      ASSUME_ALIGNED(x3);  ASSUME_ALIGNED(x4); 
 
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif       
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
    }
    double rpp2(double * Restrict y,double a1,double * Restrict x1,
		double a2,double * Restrict x2,
		double a3,double * Restrict x3,double a4,double * Restrict x4,
		double a5,double * Restrict x5,double atol,double rtol)
    {
      double ret=0;
      ASSUME_ALIGNED(x1);  ASSUME_ALIGNED(x2);
      ASSUME_ALIGNED(x3);  ASSUME_ALIGNED(x4); 
      ASSUME_ALIGNED(x5);  ASSUME_ALIGNED(y);
#ifdef ROCK4_OMP
#pragma omp parallel for
#endif       
      for(int i=0;i<size;i++)
	{
	  double ab=ABS(y[i])*rtol+atol;
	  double d=1.e0/ab;
	  double s=a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i]+a5*x5[i];
	  ret+=s*s*d*d;
	}
      return ret;
    }
 
  };
};

