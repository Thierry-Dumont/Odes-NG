#ifndef Rock4__h
#define Rock4__h
#include "Rock4Coeffs.hpp"
#include "ToolSequential.hpp"
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include "history.hpp"
#include "OdesException.hpp"
#include <cmath>
#include <string>
#include <cstdlib>
using namespace std;
///////////////////////////////////////////////////////////////////////////
/// A replacement for A. Abdulle rock4.f code for Rock4 method.
///
/// Template arguments:
///      - Fonc: Object function which defines the problem.
///      - Engine: defines sequential exectution or parallelism
///                based on Open-MP or Tbb.
//////////////////////////////////////////////////////////////////////////
namespace odes {
  template<class Fonc,class Engine=ToolSequential<Fonc> > class Rock4
  {
#ifdef Rock4_history
    history H;
#endif
    Fonc* Fp;
    const int size;
    double atol,rtol;//tolerances for time steps adaptation.
    double uround,told;
    double hlast;
    double *fn,*k1,*k2,*k3,*k4,*temp,*g0,*Vtemp;
    int mp[3],mdeg,funccal,nrej,nsteps,nacc,nstagesmax;
    Rock4Coeffs Coeffs;
 
    //! switch 2 adresses.
    //! \param a
    //! \param b
    inline void switchAdr(double *&a,double *&b)
    {
      double *c=a; a=b;b=c;
    }
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
    //! The heart of the computation; perform a dt step.
    //! \param y IN/OUT.
    //! \dt time step.
    //! \note return error estimator value.
     double kernel(double * Restrict &y,double dt)
    {
      double *g1=temp;
      Fonc& F=*Fp;// told=0.;
      double mz=mp[1],mr=mp[2];

      Engine E(mr,mz,size,dt);//E(mz,mr...)
      
      //
      //3 terms recurrence formula:
      //

      E.cl(g1,y,fn); 
      if(mdeg<2) switchAdr(g1,g0);//result in g0. yjm1 dans fortran.
      // otherwise result remains in g1==yjm1
      for(int i=2;i<=mdeg;i++)
	{
	  if(i==2)
	    {
	      E.apply(F,g0,g1,y);//g0<-h*mu F(g1)-nu*g1+kappa*y
	      switchAdr(g1,g0);
	    }
	  else
	    {
	      E.apply(F,i,g1,g0,k1);
	      switchAdr(g1,g0);
	    }
	}
      if(mdeg>=2) switchAdr(g1,g0);
      temp=g1;
      funccal+=mdeg-1;
      // result of the recurrence is always in g0.
      //
      // The second RK formula (3 steps).
      //
      //1rst step
      F(g0,k1);
 
      //2nd step:
      double teta1=dt*Coeffs.fpa(mz,1);
      E.cl(temp,g0,teta1,k1);
      F(temp,k2);
      //3rd step
      teta1=dt*Coeffs.fpa(mz,2);
      double teta2=dt*Coeffs.fpa(mz,3);
      E.cl(temp,g0,teta1,k1,teta2,k2);
      F(temp,k3);
      
      //4th step:
      teta1=dt*Coeffs.fpa(mz,4);teta2=dt*Coeffs.fpa(mz,5);
      double teta3=dt*Coeffs.fpa(mz,6);
      E.cl(temp,g0,teta1,k1,teta2,k2,teta3,k3);
      F(temp,k4);

      funccal+=4;
      //! finish:
      teta1=dt*Coeffs.fpb(mz,1);teta2=dt*Coeffs.fpb(mz,2);
      teta3=dt*Coeffs.fpb(mz,3);
      double teta4=dt*Coeffs.fpb(mz,4);
      E.ct(temp,g0,teta1,k1,teta2,k2,teta3,k3,teta4,k4);
 
      //result is temporary in temp. If the step is accepted, we need to
      // switch temp and y.
      // error estimator:
      teta1=dt*Coeffs.fpbe(mz,1)-teta1;teta2=dt*Coeffs.fpbe(mz,2)-teta2;
      teta3=dt*Coeffs.fpbe(mz,3)-teta3;teta4=dt*Coeffs.fpbe(mz,4)-teta4;
      double teta5=dt*Coeffs.fpbe(mz,5);
      Vtemp=g0;
      F(temp,Vtemp);
      
      //if the step is successfull, we must switch g0 and fn.
      ++funccal;
      double err=sqrt(E.rpp2(temp,teta1,k1,teta2,k2,teta3,k3,
			     teta4,k4,teta5,Vtemp,
			     atol,rtol)/size);
       
      //cout<<"err= "<<err<<endl;
      return err;
    }
  public:
    //! constructor
    //! \param  F RHS.
    Rock4(Fonc& F): Fp(&F),size(Fp->size())
    {
      k1=allocDoubleArray(size);k2=allocDoubleArray(size);
      k3=allocDoubleArray(size);k4=allocDoubleArray(size);
      g0=allocDoubleArray(size);fn=allocDoubleArray(size);
      temp=allocDoubleArray(size); 
      // do not allocate Vtemp.
      uround=1.e-16; atol=1.e-5; rtol=1.e-5;
      
    }
    //! destructor
    ~Rock4()
    {
      destroyDoubleArray(k1); destroyDoubleArray(k2); 
      destroyDoubleArray(k3); destroyDoubleArray(k4);
      destroyDoubleArray(g0); destroyDoubleArray(fn);
      destroyDoubleArray(temp);
    }
    //! fix values to tolerances (set by default in constructor).
    //! \param _atol
    //! \param _rtol
    void setTolerances(double _atol, double _rtol)
    {
      atol=_atol; rtol=_rtol;
    }
    //! number of stages used (we count both methods).
    inline int NbStages() const {return nstagesmax;}
    //! number of RHS computed
    inline int NbRhsComputed() const {return funccal;}
    //! number of steps performed.
    inline int NbSteps() const {return nsteps;}
    //! number of accepted steps.
    inline int NbAccepted() const {return nacc;}
    //! number of rejected steps.
    inline int NbRejected() const {return nrej;}
    //! last accepted time step.
    inline double LastAcceptedTimeStep() const {return hlast;}
    //! number of unknowns.
    inline int nbUnkn() const {return size;}
#ifdef Rock4_history
    inline history& get_history(){return H;}
#endif
    //! operator() performs the step, by delegating to kernel.
    //! \param  y IN/OUT. Initial values on input, final values on output.
    //! \param t0 initial value of time
    //! \param tend final time.
    //! \param _h proposed time step.
    //! \note We integrate from t0 to tend
    //! \note *BEWARE* adress y can be changed!
    void operator()(double *&y,double t0,double tend, double _h)
    {
      int mdego=0; nsteps=0;
      nrej=0; nacc=0;nstagesmax=0;
      double eigmax=0.0;
      double hnew; double facmax=5.0,hp=_h; 
      //bool last=false; 
      int nrho=0; int nrejloc=0;
      Fonc& F=*Fp;
      int pass=0;
      F(y,fn);
      funccal=1;
      double t=t0,h=_h,errp=-1.e0;
      told=0.0;
      bool reject=false;
 
      while(t<tend)
	{
	  if(1.01e0*h>=ABS(tend-t))
	    h=ABS(tend-t);

	  if (h<10.0e0*uround)
	    throw OdesException("Rock4: step too small h=",h,"limit:",
				     10.0e0*uround);
	  if(nrho==0)
	    eigmax=F.rho();
	  
	  mdeg=sqrt((3.e0+h*eigmax)/0.353e0)+1;

	  ++pass;

	  if (mdeg>152)
	    {
	      h=0.8e0*(152.e0*152.e0*0.353e0-3.e0)/eigmax;
	      mdeg=152;
	      //last=false;
	    }
	  mdeg=max(mdeg,5)-4;
	  if(mdeg!=mdego)
	    {
	      mdegre(mdeg);
	      mdego=mdeg;
	      if(mdeg+4>nstagesmax) nstagesmax=mdeg+4;
	    }
	  double err=kernel(y,h);

	  ++nsteps;
	  // error control:
	  double fac=sqrt(sqrt(1.e0/err));
	  if (errp>0.0 && !reject) 
	    {
	      double facp=sqrt(sqrt(errp))*fac*fac*(h/hp);
	      fac=min(fac,facp);
	    }
	  if (reject) facmax=1.e0;
	  fac=min(facmax,max(0.1e0,0.8e0*fac));
	  hnew=h*fac;

	  if(err<1.0e0)
	    {
	      //accepted step.
	      ++nacc; facmax=2.0e0;
	      t+=h;
	      hlast=h;
	      if (reject)
		{
		  hnew=min(hnew,h);
		  if (tend<t)  hnew=max(hnew,h);
		  reject=false;
		  nrejloc=0;
		}
#ifdef Rock4_history
	      H.push(h,true);
#endif
	      hp=h; h=hnew;
	      errp=err;
	      switchAdr(temp,y); switchAdr(g0,fn);
	      ++nrho;
	      nrho=(nrho+1)%10;
	      if(1.1*h>=ABS(tend-t)) h=ABS(tend-t);
	    }
	  else
	    {
#ifdef Rock4_history
	      H.push(h,false);
#endif
	      reject= true;
	      //last=false;
	      h= 0.8e0*hnew;
	      ++nrej;
	      if(nsteps==0) h/=10.e0;
	      if (told==t) 
		{
		  ++nrejloc;
		  if (nrejloc==10) h=1.0e-5;
		}
	      told=t;
	     
	      nrho= nrho!=0? 0:1;
	    }
	}
    }
   
  };
};
#endif
