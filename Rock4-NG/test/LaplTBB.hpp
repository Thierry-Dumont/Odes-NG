#ifndef Lapl1OMP__h
#define Lapl1OMP__h
#include "MacrosForCompilers.hpp"
#include <iostream>
using namespace std;
#define MIN(x,y) ((x)<(y)?(x):(y))
class ApplyLaplTBB
{
  static const int bsize=BLOCK_FACTOR/sizeof(double);
  double *x, *y;
  const int size; const double uh2;
public:
  ApplyLaplTBB(double * Restrict _x,double * Restrict  _y,int _size,double _uh2) 
    :x(_x),y(_y),size(_size),uh2(_uh2)
  {

  }
  //! Copy contructor.
  ApplyLaplTBB(const ApplyLaplTBB& A):x(A.x),y(A.y),size(A.size),uh2(A.uh2)
  {
  }
   void operator()(const blocked_range<size_t>& r) const
    {
      int i1=r.begin()*bsize,i2=MIN(size,r.end()*bsize);
      bool bd=i2==size;
      if(i1==0)
	{
	  y[0]= uh2*(x[1]-x[0]); ++i1;
	}
      if(bd) --i2;
      for(int i=i1;i<i2;i++)
	y[i]=uh2*(x[i-1]-2.0*x[i]+x[i+1]);
      if(bd)
	y[size-1]= uh2*(x[size-2]-x[size-1]);
    }
};
class LaplTBB
{
  static const int bsize=BLOCK_FACTOR/sizeof(double);
  const int Size;
  double nu,h,uh2,rspec;
  int sizep;
public:
  LaplTBB(int _Size): Size(_Size)
  {
    
    nu=0.01;// diffusion.
    h=1./(Size-1.e0); uh2=nu/(h*h); rspec=2.*uh2;
    sizep=Size/bsize+ (Size%bsize? 1:0);
    cout<<"rspec= "<<rspec<<endl;
  }
  ~LaplTBB(){}
 
  inline void operator()(double * Restrict x,double * Restrict  y) const
  {
    ASSUME_ALIGNED(x); ASSUME_ALIGNED(y);
    parallel_for(blocked_range<size_t>(0,sizep),
		 ApplyLaplTBB(x,y,Size,uh2));
  }
  inline  double rho() const {
    return rspec;
  } 
  inline int size() const {return Size;}
};

#endif
