#ifndef ToolTbb__h
#define ToolTbb__h
#include "Rock4Coeffs.hpp"
#include "OdesException.hpp"
#include "MacrosForCompilers.hpp"
#include "TbbUtils.hpp"
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
#define P2(x)  ((x*x))
#include <iostream>
using namespace std;
/////////////////////////////////////////////////////////////////////////////
/// Perform basic operations for Rock4. 
/// Parallelism using tbb.
/////////////////////////////////////////////////////////////////////////////
namespace odes {
  template<class Fonc> class ToolTbb
  {
    static const int bsize=BLOCK_FACTOR/sizeof(double);
    int mr,mz,size,sizep; double h;
    Rock4Coeffs Coeffs;
  public:
    void print(double x[],string s="")
    {
      std::cout.precision(16);
      cout<<endl<<s<<" ";
      for(int i=0;i<size;i++)
	cout<<scientific<<x[i]<<" ";cout<<endl;
    }
    //! Constructor
    //! \param _mr see Rock4.hpp
    //! \param _mz see Rock4.hpp
    //! \param _size see Rock4.hpp
    //! \param h time step.
    ToolTbb(int _mr,int _mz,int _size,double _h):
      mr(_mr),mz(_mz),size(_size),h(_h)
    {
      sizep=size/bsize+ (size%bsize? 1:0);
    }
    //! \Destructor
    ~ToolTbb(){}
    void cl(double * Restrict g,double * Restrict y,double * Restrict fn)
    {
      double temp=h*Coeffs.recf(mr);
      parallel_for(blocked_range<size_t>(0,sizep),
		   DoLinearComb(size,g,y,temp,fn));
    }
    void cl(double * Restrict cible,double * Restrict y,double a1,double x1[])
    {
      parallel_for(blocked_range<size_t>(0,sizep),
		   DoLinearComb(size,cible,y,a1,x1));
    }
    void cl(double  * Restrict cible,double  * Restrict y,double a1,
	    double  * Restrict x1,double a2,double  * Restrict x2)
    {
      parallel_for(blocked_range<size_t>(0,sizep),
		   DoLinearComb(size,cible,y,a1,x1,a2,x2));
    }
    void cl(double * Restrict cible,double * Restrict y,double a1,
	    double * Restrict x1,double a2,double * Restrict x2,
	    double a3,double * Restrict x3)
    {
      parallel_for(blocked_range<size_t>(0,sizep),
		   DoLinearComb(size,cible,y,a1,x1,a2,x2,a3,x3)
		   );
    }
    void ct(double cible[],double y[],double a1,double x1[],
	    double a2,double x2[],double a3,double x3[],double a4,double x4[])
    {
      parallel_for(blocked_range<size_t>(0,sizep),
		   DoLinearComb(size,cible,y,a1,x1,a2,x2,a3,x3,a4,x4)
		   );
    }
    double rpp2(double * Restrict y,double a1,double * Restrict x1,
		double a2,double * Restrict x2,
		double a3,double * Restrict x3,double a4,double * Restrict x4,
		double a5,double * Restrict x5,double atol,double rtol)
    {
      RPP2 R(size,y,a1,x1,a2,x2,a3,x3,a4,x4,a5,x5,atol,rtol);
      parallel_reduce(blocked_range<size_t>(0,sizep),R);
      return R.result;

    }
    void apply(Fonc& F,int deg,double g1[], double g0[], double vtemp[])
    {
      double temp1=h*Coeffs.recf(mr+2*(deg-2)+1),
	temp3=-Coeffs.recf(mr+2*(deg-2)+2);
      double temp2=1.0-temp3;
      F(g1,vtemp);
      // parallel_for(blocked_range<size_t>(0,size),
      // 		   Apply<Fonc>(F,temp1,temp2,temp3,g0,g1)
      // 		   );
      parallel_for(blocked_range<size_t>(0,sizep),
      		   DoLinearComb(size,g0,temp1,vtemp,temp2,g1,temp3,g0)
      		   );
    }
    void apply(Fonc& F,double g2[], double g1[],double g0[])
    {
      // here, degree==2.
      double temp1=h*Coeffs.recf(mr+1),	temp3=-Coeffs.recf(mr+2);
      double temp2=1.0-temp3;
      F(g1,g2);
      parallel_for(blocked_range<size_t>(0,sizep),
      		   DoLinearComb(size,g2,temp1,g2,temp2,g1,temp3,g0)
      		   );  
    }
  };
};
#endif
