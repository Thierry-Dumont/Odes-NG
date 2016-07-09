#ifndef Rock4__h
#define Rock4__h
#include "Rock4Coeffs.hpp"
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include "GenericException.hpp"
#include <cmath>
#include <cstdlib>
using namespace std;
namespace odes{
///////////////////////////////////////////////////////////////////////////
/// A replacement for A. Abdulle rock4.f code for Rock4 method.
///  *We treat here linear and affine case* ie RHS must be:
///        x-> Ax
///    or 
///        x-> Ax+ B with B independent of time.
///
/// Template arguments:
///      - Fonc: Object function which defines the problem.
///      - Affine: true iff we treat the affine case.
///              
///  The class Fonc must define the following public methods:
///      1) contructor:  
///                Fonc(int n) n: size (number of unknowns).
///      2) destructor.
///      3)  double rho()
///            return spectral radius of A: 
///      4) int size()        
///            return n, the size of vectors (see contructor).
///      5) void operator()(double * in, double * s, 
///                            double a, double b,double *out)
///           compute out[]= A*a* in[] +b s[]
///      6) void operator()(double * in, double * s, 
///                             double a, double b, double *g,double *out)
///           compute out[]=A*a*in[]+ b*s[] +g[]
///      7) void operator()(double * in, double * s,
///                              double a, double b,double c,double *out)
///           compute out[]==A*a*in[]+ b*in[] + c*s[]
///      8)  void operator()(double * in, double * s, 
///                              double a, double b,double c,double *g,
///                              double *out)
///           compute out[]==A*a*in[]+ b*in[] + c*s[]+ g[]
///
///    where in[], s[], g[] and out[] are vectors of size n.        
//////////////////////////////////////////////////////////////////////////
namespace odes {
  template<class Fonc,bool Affine=false> class Rock4L
  {
    Fonc* Fp;
    const int size;
    double atol,rtol;//tolerances for time steps adaptation.
    double uround,told;
    double t; 
    double *fn,*k1,*k2,*temp,*g0,*g1,*g2;
    double *BigArray;
    int mp[3],mdeg,funccal,nrej,nsteps,nacc,nstagesmax,sizesplit;
    Rock4Coeffs Coeffs;
    double ac[5][4],bc[5],hc[5];

    //! We want to compute polynomial for the 2nd formula in Horner form.
    //! so we must compute the coefficients in this form.
    void makeHornerCoeffs()
    {
      hc[4]=ac[2][1]*ac[3][2]*ac[4][3]*bc[4]; 
      hc[3]=ac[2][1]*ac[3][2]*bc[3] + ac[2][1]*ac[4][2]*bc[4] + 
	ac[3][1]*ac[4][3]*bc[4] + ac[3][2]*ac[4][3]*bc[4];
      hc[2]= ac[2][1]*bc[2] + ac[3][1]*bc[3] + ac[3][2]*bc[3] + 
	ac[4][1]*bc[4] + ac[4][2]*bc[4] + ac[4][3]*bc[4];
      hc[1]= bc[1] + bc[2] + bc[3] + bc[4];
      hc[0]=1.;
    } 
    //! make sub arrays.
    inline void splitBigArray()
    {
      k1=BigArray; g0=k1+size; g1=g0+size; g2=g1+size;
    }
    //! permutate 3 adresses g0,g1,g2.
    //! \param g0 
    //! \param g1 
    //! \param g2
    inline void permut3(double *&g0,double *&g1,double *&g2)
    {
      double *temp=g0;
      g0=g1; g1=g2;g2=temp;
    }
    //! switch 2 adresses.
    //! \param a
    //! \param b
    inline void switchAdr(double *&a,double *&b)
    {
      double *c=a; a=b;b=c;
    }
    //! find degree of formula.
    void mdegre(int& mdeg)
    {
      mp[2]=1;
      //-------- Find the degree.--------      
      for(int i=1;i<=50;i++)
	{
	  if ((Coeffs.ms(i)/mdeg)>=1)
	    {
	      mdeg=Coeffs.ms(i);
	      mp[1]=i;
	      break;
	    }
	  else
	    mp[2]=mp[2]+Coeffs.ms(i)*2-1;
	}
    }
    //! The heart of the computation; perform a dt step, in the linear
    //! case.
    //! \param y IN/OUT.
    //! \param dt time step.
    //! \note we solve dy/dt= Ay
     void kernel1(double *y,double dt)
    {
      Fonc& F=*Fp;// told=0.;
      double mz=mp[1],mr=mp[2];

      splitBigArray();
      //
      //3 terms recurrence formula:
      //
      F(y,y,dt*Coeffs.recf(mr) ,1.0,g1);
 
      for(int i=2;i<=mdeg;i++)
	{
	  if(i==2)
	    {
	      double temp1=dt*Coeffs.recf(mr+1),temp3=-Coeffs.recf(mr+2);
	      double temp2=1.0-temp3;
	      F(g1,y,temp1,temp2,temp3,g2);
	      permut3(g0,g1,g2);

	    }
	  else
	    {
	      double temp1=dt*Coeffs.recf(mr+2*(i-2)+1);
	      double temp3=-Coeffs.recf(mr+2*(i-2)+2);
	      double temp2=1.0-temp3;
	      F(g1,g0,temp1,temp2,temp3,g2);
	      permut3(g0,g1,g2);
	    }
	}
      switchAdr(g1,g0);
      temp=g1;
      funccal+=mdeg+1;
      ++nsteps;
      //
      // result of the recurrence is always in g0.
      //
      // The second RK formula (4 steps).
      //
      //coefficients:
      
      ac[2][1]=dt*Coeffs.fpa(mz,1);
      ac[3][1]=dt*Coeffs.fpa(mz,2);ac[3][2]=dt*Coeffs.fpa(mz,3);
      ac[4][1]=dt*Coeffs.fpa(mz,4);ac[4][2]=dt*Coeffs.fpa(mz,5);
      ac[4][3]=dt*Coeffs.fpa(mz,6);
      bc[1]=dt*Coeffs.fpb(mz,1);bc[2]=dt*Coeffs.fpb(mz,2);
      bc[3]=dt*Coeffs.fpb(mz,3);bc[4]=dt*Coeffs.fpb(mz,4);
      makeHornerCoeffs();

      //
 
      F(g0,g0,hc[4],hc[3],k1);
      F(k1,g0,1.0,hc[2],temp);
      F(temp,g0,1.0,hc[1],k1);
      F(k1,g0,1.0,hc[0],y);
      //
      funccal+=4;

    }
    //! The heart of the computation; perform a dt step, in the affine
    //! case.
    //! \param y IN/OUT.
    //! \param b the constant in time part of the RHS
    //! \param dt time step.
    //! \note we solve dy/dt= Ay + B, with a constant B.
    void kernel2(double *y,double *b,double dt)
    {
   
      //g1=temp;
      Fonc& F=*Fp;// told=0.;
      double mz=mp[1],mr=mp[2];

      splitBigArray();
      //
      //3 terms recurrence formula:
      //
      F(y,y,dt*Coeffs.recf(mr) ,1.0,b,g1);
 
      for(int i=2;i<=mdeg;i++)
	{
	  if(i==2)
	    {
	      double temp1=dt*Coeffs.recf(mr+1),temp3=-Coeffs.recf(mr+2);
	      double temp2=1.0-temp3;
	      F(g1,y,temp1,temp2,temp3,b,g2);
	      permut3(g0,g1,g2);

	    }
	  else
	    {
	      double temp1=dt*Coeffs.recf(mr+2*(i-2)+1);
	      double temp3=-Coeffs.recf(mr+2*(i-2)+2);
	      double temp2=1.0-temp3;
	      F(g1,g0,temp1,temp2,temp3,b,g2);
	      permut3(g0,g1,g2);
	    }
	}
      switchAdr(g1,g0);
      temp=g1;
      funccal+=mdeg+1;
      ++nsteps;
      //
      // result of the recurrence is always in g0.
      //
      // The second RK formula (4 steps).
      //
      //coefficients:
      
      ac[2][1]=Coeffs.fpa(mz,1);
      ac[3][1]=Coeffs.fpa(mz,2);ac[3][2]=Coeffs.fpa(mz,3);
      ac[4][1]=Coeffs.fpa(mz,4);ac[4][2]=Coeffs.fpa(mz,5);
      ac[4][3]=Coeffs.fpa(mz,6);
      bc[1]=Coeffs.fpb(mz,1);bc[2]=Coeffs.fpb(mz,2);
      bc[3]=Coeffs.fpb(mz,3);bc[4]=Coeffs.fpb(mz,4);
      makeHornerCoeffs();
      double dt2=dt*dt, dt3= dt2*dt;
      //
      F(g0,b,dt,dt,k1);// dt A g0+ dt b -> k1

      // now we have a polynomial of degree 3 to evaluate: 
      F(k1,k1,hc[4]*dt3,hc[3]*dt2,g1);
      F(g1,k1,1.0,hc[2]*dt,g2);
      F(g2,k1,1.0,hc[1],g1);
      // finish:
      ASSUME_ALIGNED(y);ASSUME_ALIGNED(g0);ASSUME_ALIGNED(g1);
#include "Ivdep.hpp"
      for(int i=0;i<size;i++)
	y[i]=g0[i]+g1[i];
      
      funccal+=4;

    }

  public:
    //! constructor
    //! \param  F RHS.
    Rock4L(Fonc& F): Fp(&F),size(Fp->size())
    {
      int bk=BLOCK_FACTOR/sizeof(double);
      sizesplit=BLOCK_FACTOR*(size/bk+ size%bk==0? 0:1);
      BigArray=allocDoubleArray(4*size);//Allocate(4*size);
      splitBigArray();
       // do not allocate Vtemp. 
      uround=1.e-16;
    }
    //! destructor
    ~Rock4L()
    {
      //Suppress(BigArray);
      destroyDoubleArray(BigArray);
    }
    
    //! number of stages used (we count both methods).
    inline int NbStages() const {return nstagesmax;}
    //! number of RHS computed
    inline int NbRhsComputed() const {return funccal;}
    //! number of steps performed.
    inline int NbSteps() const {return nsteps;}
    //! return a correctly aligned vector, for unknowns.
    //inline double* alloc_vec() {return Allocate(size);} 
    //! delete a vector allocated with alloc_vec()
    //inline void delete_vec(double * v){Suppress(v);}
    //! number of unknowns.
    inline int nbUnkn() const {return size;}
    //! value of "time":
    inline double Time() const {return t;}
    
    //! operator() performs the step IN THE LINEAR CASE, delegating to kernel1.
    //! \param  y IN/OUT. Initial values on input, final values on output.
    //! \param t0 initial value of time
    //! \param tend final time.
    //! \param _h proposed time step.
    //! \note We integrate from t0 to tend
    //! \note  we solve dy/dt= Ay.
    void operator()(double *y,double t0,double tend, double _h)
    {
      if(Affine)
	throw GenericException("Rock4::operator() called as affine");
      int mdego=0; nsteps=0;
      nrej=0; nacc=0;nstagesmax=0;
      Fonc& F=*Fp;
      int pass=0;
      funccal=0;
      t=t0;
      double h=_h;
      told=0.0;
      double eigmax=F.rho();
      
      while(t<tend)
	{
	  if(1.01e0*h>=abs(tend-t))
	    h=abs(tend-t);

	  if (h<10.0e0*uround)
	    throw GenericException("Rock4: step too small h=",h,"limit:",
				     10.0e0*uround);

	  ++pass;
	  mdeg=sqrt((3.e0+h*eigmax)/0.353e0)+1;
	  if (mdeg>152)
	    throw GenericException("Rock4: degree limited to 152."
				   "Here you need:",mdeg,
				   "You should decrease time step");
	  mdeg=max(mdeg,5)-4;
	  if(mdeg!=mdego)
	    {
	      mdegre(mdeg);
	      mdego=mdeg;
	      if(mdeg+4>nstagesmax) nstagesmax=mdeg+4;
	    }
	  kernel1(y,h);
	  t+=h;
	  if(t+h>tend) h=tend-t;
	}
    }
    //! operator() performs the step IN THE AFFINE CASE, delegating to kernel1.
    //! \param  y IN/OUT. Initial values on input, final values on output.
    //! \param b RHS.
    //! \param t0 initial value of time
    //! \param tend final time.
    //! \param _h proposed time step.
    //! \note We integrate from t0 to tend
    //! \note we solve dy/dt= Ay+b.
    void operator()(double *y,double *b, double t0,double tend, double _h)
    {
      if(!Affine)
	throw GenericException("Rock4::operator() called as not affine");
      int mdego=0; nsteps=0;
      nrej=0; nacc=0;nstagesmax=0;
      Fonc& F=*Fp;
      int pass=0;
      funccal=0;
      t=t0;
      double h=_h;
      told=0.0;
      double eigmax=F.rho();
      
      while(t<tend)
	{
	  if(1.01e0*h>=abs(tend-t))
	    h=abs(tend-t);

	  if (h<10.0e0*uround)
	    throw GenericException("Rock4: step too small h=",h,"limit:",
				     10.0e0*uround);

	  ++pass;
	  mdeg=sqrt((3.e0+h*eigmax)/0.353e0)+1;
	  if (mdeg>152)
	    throw GenericException("Rock4: degree limited to 152."
				   "Here you need:",mdeg,
				   "You should decrease time step");
	  mdeg=max(mdeg,5)-4;
	  if(mdeg!=mdego)
	    {
	      mdegre(mdeg);
	      mdego=mdeg;
	      if(mdeg+4>nstagesmax) nstagesmax=mdeg+4;
	    }
	  kernel2(y,b,h);
	  t+=h;
	  if(t+h>tend) h=tend-t;
	}
    }
  };
};
};
#endif
