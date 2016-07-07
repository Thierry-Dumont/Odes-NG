#ifndef ToolOMP__h
#define ToolOMP__h
#include "Rock4Coeffs.hpp"
#include "OdesException.hpp"
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
#define P2(x)  ((x*x))
#include <iostream>
using namespace std;
namespace odes {
  ///////////////////////////////////////////////////////////////////////////
  /// Perform basic operations for Rock4. Parallelism with OMP.
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
	cout<<scientific<<x[i]<<" ";cout<<endl;
    }
    ToolOMP(int _mr,int _mz,int _size,double _h):
      mr(_mr),mz(_mz),size(_size),h(_h){}
    ~ToolOMP(){}
    //! The recurrence for one degree.
    //! \param F RHS.
    //! \param deg the degree
    //! \param g1 preceeding term
    //! \param g0 on input: 2 terms ago; result on output.
    void apply(Fonc& F,int deg,double g1[], double g0[])
    {
      double temp1=h*Coeffs.recf(mr+2*(deg-2)+1),
	temp3=-Coeffs.recf(mr+2*(deg-2)+2);
      double temp2=1.0-temp3;

#pragma omp parallel for
      for(int i=0;i<size;i++)
	g0[i]=temp1*F(i,g1)+temp2*g1[i]+temp3*g0[i];
      
    }
    void apply(Fonc& F,int deg,double g2[], double g1[],double g0[])
    {
      // print(g0,"y dans apply");
      // print(g2,"g0 dans apply");
      // print(g1,"g1 dans apply");
      double temp1=h*Coeffs.recf(mr+2*(deg-2)+1),
	temp3=-Coeffs.recf(mr+2*(deg-2)+2);
      double temp2=1.0-temp3;
     
#pragma omp parallel for
      for(int i=0;i<size;i++)
	g2[i]=temp1*F(i,g1)+temp2*g1[i]+temp3*g0[i];
    }
    void cl(double *g,double *y,double *fn)
    {
      double temp=h*Coeffs.recf(mr);
#pragma omp parallel for
      for(int i=0;i<size;i++)
	g[i]=y[i]+temp*fn[i];
    }
    void cl(double cible[],double y[],double a1,double x1[])
    {
#pragma omp parallel for
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i];
    }
    void cl(double cible[],double y[],double a1,double x1[],
	    double a2,double x2[])
    {
#pragma omp parallel for
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i]+a2*x2[i];
    }
    void cl(double cible[],double y[],double a1,double x1[],
	    double a2,double x2[],double a3,double x3[])
    {
#pragma omp parallel for
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i];
    }
    void ct(double cible[],double y[],double a1,double x1[],double a2,double x2[],
	    double a3,double x3[],double a4,double x4[])
    {
#pragma omp parallel for
      for(int i=0;i<size;i++)
	cible[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
    }
    double rpp2(double y[],double a1,double x1[],double a2,double x2[],
		double a3,double x3[],double a4,double x4[],
		double a5,double x5[],double atol,double rtol)
    {
      double ret=0;

#pragma omp parallel for reduction(+:ret)
      for(int i=0;i<size;i++)
	ret+=P2((a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i]+a5*x5[i])/
		(ABS(y[i])*rtol+atol));
      return ret;
    }
  };
};
#endif
