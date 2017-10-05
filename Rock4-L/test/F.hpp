#include <cmath>
#include <string>
#include <iostream>
#include "MacrosForCompilers.hpp"
using namespace std;
class F
{
public :
  int Size;
  double h,eps;
  //! spectral radius
  double rho ()
  {
     double pi=3.14159265359;
     double R=0;
     for(int i = 0; i<Size; i++)
       R=max(abs(  (2/(h*h))*(cos(i*pi/(Size))-1)),R);
 
    return R ;
  }
  //! constructor
  //! \param n size (number of points, including boundaries).
  F (int n)
  {
    Size=n;
    h=1./((double)(n)-1);
   }
  //! destructor
  ~F(){}
  //! return number of points
  int size ()const
  {
    return Size;
  }



  //! a \Delta in  + b s. 
  void operator()(double * in, double * s, double a, double b,double *out)
  {
    double h1=h*h,ah1=a/h1;
    out[0]=(-2.*in[0]+2.*in[1])*ah1 +b*s[0];
    ASSUME_ALIGNED(out);ASSUME_ALIGNED(in);ASSUME_ALIGNED(s);
#pragma omp parallel for
    for(int i=1;i<Size-1;i++)
      {
	out[i]=in[i-1]-2*in[i]+in[i+1];
	out[i]=out[i]*ah1+b*s[i];
      }
    out[Size-1]=(2.*in[Size-2]-2.*in[Size-1])*ah1+b*s[Size-1];
  

  }
 //! a (\Delta in + g)  + b s 
  void operator()(double * in, double * s, double a, double b, double *g,
		  double *out)
  {
    double h1=h*h;
   
    out[0]=a*((-2.*in[0]+2.*in[1])/h1 +g[0])+b*s[0];
    ASSUME_ALIGNED(out);ASSUME_ALIGNED(in);ASSUME_ALIGNED(g);
#pragma omp parallel for
    for(int i=1;i<Size-1;i++)
    {
      out[i]=in[i-1]-2*in[i]+in[i+1];
      out[i]=a*(out[i]/h1+g[i])+b*s[i];
    }
    out[Size-1]=a*((2.*in[Size-2]-2.*in[Size-1])/h1+g[Size-1])
      +b*s[Size-1];
  }
  // out = aAin + b in + c s
  void operator()(double * in, double * s, double a, double b,double c,
		  double *out)
  {
    ASSUME_ALIGNED(out);ASSUME_ALIGNED(in);ASSUME_ALIGNED(s);
    double h1=h*h,ah1=a/h1; 
 
    out[0]=(-2.*in[0]+2.*in[1])*ah1;
#pragma omp parallel for
    for(int i=1;i<Size-1;i++)
      {
	out[i]=in[i-1]-2*in[i]+in[i+1];
	out[i]=out[i]*ah1;
      }
    out[Size-1]=(2.*in[Size-2]-2.*in[Size-1])*ah1;
    
#pragma omp parallel for
    for(int i=0;i<Size;i++)
    {  
      out[i]=out[i]+b*in[i]+c*s[i];
    }

  }
  // out = a(Ain+g) + b in + c s
  void operator()(double * in, double * s, double a, double b,double c,
		  double *g,double *out)
  {
    ASSUME_ALIGNED(out);ASSUME_ALIGNED(in);
    ASSUME_ALIGNED(s);ASSUME_ALIGNED(g);
    
    double h1=h*h; 
 
    out[0]=a*((-2.*in[0]+2.*in[1])/h1+g[0])+b*in[0]+c*s[0];

#pragma omp parallel for
    for(int i=1;i<Size-1;i++)
      {
	out[i]=a*((in[i-1]-2*in[i]+in[i+1])/h1+g[i])+b*in[i]+c*s[i];
      }
    out[Size-1]=a*((2.*in[Size-2]-2.*in[Size-1])/h1+g[Size-1])+b*in[Size-1]
    +c*s[Size-1];
  
  }

};
