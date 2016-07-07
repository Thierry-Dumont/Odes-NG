#ifndef rock4__h
#define rock4__h
#include "OdesException.hpp"
#include "Rock4Coeffs.hpp"
#include <iostream>
#include <utility>
#include <math.h>
#include "tbb/tbb.h"
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
#define P2(x)  ((x*x))
using namespace std;
using namespace tbb;
///////////////////////////////////////////////////////////////////////
///
/// Class Vect: variables.
///    *public methods of Vect: 
///        nunc(): number of variables
///        double operator()(int i): return variable i !Z! i in [1,n()]
///
/// Class Fct: second members
///     * public (or protected) methods of Fct:
///         void operator(Vect& Y, const Vect& X): compute Y=Fct(X)
///         constructor: F(Vect& X)
///         double rho(const Vect& X): estimated Jacobian spectral
///                                                    radius at X. 
///////////////////////////////////////////////////////////////////////
template<class Vect,class Fct> class rock4: private Rock4Coeffs 
{
  bool initialized;
  int n;
  Vect *Yn,*Yjm1,*Yjm2,*Yjm3,*Yjm4,*Fn, *Fnt;
  int kpass;
  int funccall;

  double eigmax;//rayon spectral
  double t;
  int mp[3];
  double hp,told,hnew,errp,err;
  double uround,rtol,atol;
  double h;
  int nrho,nrej,idid;
  bool reject,firstStepWhasRejected;
  int rejected;
  //! core of the method.
  double rfstep(Fct& F,Vect& y,int mdeg)
  {
    Vect& yn=*Yn;
    Vect& yjm1=*Yjm1; Vect& yjm2=*Yjm2;
    Vect& yjm3=*Yjm3; Vect& yjm4=*Yjm4;
    Vect& fn=*Fn; Vect& fnt=*Fnt;
    // 
    double err=0.0;
    hp=h; told=0.0;

    int mz=mp[1],mr=mp[2];
    //-----1rst stage:
    double temp1=h*recf(mr),temp2,temp3,temp4,temp5;
    double ci1=t+temp1,ci2=ci1,ci3=t;
    yjm2=yn;
 
    
    yjm1.put(1.0,yn,temp1,fn);

    
    
    if(mdeg<2)// switch!
      y=yjm1;
    //-----2nd stage:
    for(int i=2;i<=mdeg;i++)
      {
	temp1=h*recf(mr+2*(i-2)+1);
        temp3=-recf(mr+2*(i-2)+2);
        temp2=1.0-temp3;



	F(ci1,yjm1,y);
	++funccall;

	ci1=temp1+temp2*ci2+temp3*ci3;
	y.put(temp1,y,temp2,yjm1,temp3,yjm2);
	if(i<mdeg)
	  {
	    yjm2=yjm1;
	    yjm1=y;
	  }
	ci3=ci2; ci2=ci1;
      }
    
    //The finishing procedure (4-stage method).
    //Stage 1.
    temp1=h*fpa(mz,1);
    F(ci1,y,yjm1); ++funccall;

    yjm3.put(1.0,y,temp1,yjm1);
    
    //Stage 2.
    ci2=ci1+temp1;    temp1=h*fpa(mz,2);    temp2=h*fpa(mz,3);
    F(ci2,yjm3,yjm2); ++funccall;
    yjm4.put(1.0,y,temp1,yjm1,temp2,yjm2);
    //Stage 3.
    ci2=ci1+temp1+temp2;
    temp1=h*fpa(mz,4);    temp2=h*fpa(mz,5);    temp3=h*fpa(mz,6);
    F(ci2,yjm4,yjm3); ++funccall;
    fnt.put(1.0,y,temp1,yjm1,temp2,yjm2,temp3,yjm3);
    //Stage 4.
    ci2=ci1+temp1+temp2+temp3;
    temp1=h*fpb(mz,1);    temp2=h*fpb(mz,2);    temp3=h*fpb(mz,3);
    temp4=h*fpb(mz,4);
    F(ci2,fnt,yjm4); ++funccall;
    y.put(1.0,y,temp1,yjm1,temp2,yjm2,temp3,yjm3,temp4,yjm4);
    
    //Error evaluation (embedded method of order 3)
    temp1=h*fpbe(mz,1)-temp1;    temp2=h*fpbe(mz,2)-temp2;
    temp3=h*fpbe(mz,3)-temp3;    temp4=h*fpbe(mz,4)-temp4;
    temp5=h*fpbe(mz,5);
    F(t+h,y,fnt); ++funccall;
    //Atol and rtol are scalar.
    
    err=0.0;

    err=y.RPP2(temp1,yjm1,temp2,yjm2,temp3,yjm3,temp4,yjm4,temp5,fnt,
	       atol,rtol);
    err=sqrt(err/n);
 
    return err;
  }
  void mdegre(int mdeg)
  {
    mp[2]=1;
    //-------- Find the degree.--------      
    for(int i=1;i<=50;i++)
      {
	if ((ms(i)/mdeg)>=1)
	  {
	    mdeg=ms(i);
	    mp[1]=i;
	    break;
	  }
	
	mp[2]=mp[2]+ms(i)*2-1;
      }
  }
public:
  //! constructor
  rock4(double _rtol,double _atol):Rock4Coeffs(),rtol(_rtol),atol(_atol)
  {
    initialized=false;
  }
  //! constructor
  rock4(Vect& y,double& _t, double _rtol,double _atol):Rock4Coeffs(),
							      rtol(_rtol),
							      atol(_atol)
  {
    n=y.nunc();
    Yn=new Vect(y); 
    Yjm1=new Vect(y); Yjm2=new Vect(y);
    Yjm3=new Vect(y); Yjm4=new Vect(y);
    Fn=new Vect(y); Fnt=new Vect(y);
    kpass=0; funccall=0;
    t=_t;

    uround=1.e-16;
    initialized=true;
  }
  //! initialize:
  void init(Vect& y,double& _t)
  {
    n=y.nunc();
    Yn=new Vect(y,1); 
    Yjm1=new Vect(y,2); Yjm2=new Vect(y,3);
    Yjm3=new Vect(y,4); Yjm4=new Vect(y,5);
    Fn=new Vect(y,6); Fnt=new Vect(y,7);
    
    kpass=0; funccall=0;
    t=_t;

    uround=1.e-16;
    initialized=true;
  }
  //! initialized?
  inline bool is_initialized() const {return initialized;}
  //! destructor
  ~rock4()
  {
    if(initialized)
      {
	delete Yn; 
	delete Yjm1;   delete Yjm2;     delete Yjm3;     delete Yjm4; 
	delete Fn;     delete Fnt;
      }
  }
  void operator()(Fct& F,Vect& y, double tend,double _h,double _t=0.0)
  {
    if(!initialized)
      throw OdesException("rock4:: operator(): not initialized");
    kpass=0; funccall=0;
    firstStepWhasRejected=false;
    rejected=0;
    Vect& fn=*Fn; Vect& fnt=*Fnt; Vect& yn=*Yn;
    h=_h; t=_t;
    bool last=false; // bool arret=false;
    int mdego=0;
    int nsteps=0;
    double facmax=5.0;
    nrho=0;
    hp=h;
    last=false;
    reject=false; 
    err=0.0;
    nrho=0;
    idid=1;
    while(!last)
      {
	if(loglevel>2)
	  cout<<endl<<"Rock4: t= "<<t<<" tend= "<<tend<<" h= "<<_h<<endl;
	yn=y;
	if(kpass++)
	  fn=fnt;
	else
	  {
	    F(t,yn,fn);
	    ++funccall;
	  }

	bool reject=true;
	
	while(reject)
	  {
	    if(1.01*h>=ABS(tend-t))
	      {
		h=ABS(tend-t);
		last=true;
	      }
	    if (h<10.0*uround)
	      throw OdesException("rock4: step too small h=",h,"limit:",
				     10.0*uround);

	    eigmax=F.rho(y.whichIsUsed());
	    int mdeg=sqrt((3.+h*eigmax)/0.353)+1;
	    if(loglevel>2)
	      cout<<"eigmax= "<<eigmax<<" mdeg "<<mdeg<<endl;
	    if (mdeg>152)
	      {
		h=0.8*(152.*152.*0.353-3.)/eigmax;
		mdeg=152;
		last=false;
	      }
	    mdeg=max(mdeg,5)-4;
	    if (mdeg!=mdego) mdegre(mdeg);
	    //if(loglevel>0) cout<<"mdegre: "<<mdeg<<endl;
	    //Computation of an integration step.--------
	    double err=rfstep(F,y,mdeg);


	    mdego=mdeg;
	    //Error control procedure.--------
	    double fac=sqrt(sqrt((1./err)));
	    if (errp!=0. && !reject) 
	      {
		double facp=sqrt(sqrt(errp))*fac*fac*(h/hp);
		fac=min(fac,facp);
	      }
	    if (reject) facmax=1.;
	    fac=min(facmax,max(0.1,0.8*fac));
	    hnew=h*fac;
	    //Accepted step.--------

	    if(err<1.0)
	      {
		facmax=2.0;
		t+=h;
		if (reject)
		  {
		    hnew=min(hnew,h);
		    if (tend<t)  hnew=max(hnew,h);
		    reject=false;
		    if(loglevel>2)
		      cout<<"step accepted: "<<hnew<<endl;
		    //char az; cin>> az;
		    nrej=0;
		  }
		reject=false;
		hp=h;
		h=hnew;
		nrho=nrho+1; nrho=(nrho+1)%10;

	      }
	    else
	      {
		if(loglevel>2)
		  cout<<"rock4: step rejected"<<endl;
		reject= true;
		last=false;
		h= 0.8*hnew;
		++rejected;
		if(nsteps==0) h/=10.;
		if (told==t) 
		  {
		    nrej=nrej+1;
		    if (nrej==10) h=1.0e-5;
		  }
		told=t;
	      }
	  }
      }// end while here.
    

  }
  //! return spectral radius.
  inline double sp_rad() const {return eigmax;}
  //! return last time step used.
  inline double last_dt() const {return h;}
  //! return number of F evaluations.
  inline int FCalled() const {return funccall;}
};
#endif
