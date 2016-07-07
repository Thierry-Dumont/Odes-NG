#ifndef RADAU5CC__H
#define RADAU5CC__H
#include "OdesException.hpp"
#include "fortranArray.hpp"
#include "fortranVector.hpp"
#include "protos_lapack.hpp"
#include "Matrices.hpp"
#include "compat.hpp"
#include <cmath>
#include <utility>
//#include <type_traits>
#ifdef LOGRADAU5
#include "logger.hpp"
#endif
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
#define SIGN(a) (((a) >= (0.0)) ? (1.) : (-1.))
namespace odes
{
  /////////////////////////////////////////////////////////////////////
  /// Class Radau5cc
  ///
  /// A C++ implementation of (most part of) Hairer & Wanner Radau5.
  /// program.
  /// \brief A C++ implementation of Radau5.
  ////////////////////////////////////////////////////////////////////
  template<class Fonct> class Radau5cc: 
    private Matrices<(Fonct::n-Fonct::nsub)==1&&(Fonct::n-Fonct::nsup)==1,
      Fonct::Hessenberg,Fonct::n,Fonct::nsub,Fonct::nsup>,
      compat<(Fonct::n-Fonct::nsub)==1&&(Fonct::n-Fonct::nsup)==1,
      Fonct::Hessenberg>
  
  {
    static const int n=Fonct::n;
    static const int ninf=Fonct::nsub;
    static const int nsup=Fonct::nsup;
    static const bool MatrixFull=(n-ninf)==1&&(n-nsup)==1;
    static const bool Hessenberg=Fonct::Hessenberg;
    static const bool ComputeJacobianNumerically=
      Fonct::ComputeJacobianNumerically;
  public:
    typedef typename Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::MatrixReal
    MatrixReal;
  private:
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::Ibegin;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::Iend;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::calhes;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::decomr;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::decomc;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::solvereal;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::solvecomplex;
#ifdef LOGRADAU5
    logger Logger;
#endif
    // We try to keep most of the notation of Hairer & Wanner.
    enum OutNewtonLoop { converged, willNotConverge,didNotConverge };
    double SQ6,C1,C2,C1M1,C2M1,C1MC2,DD1,DD2,DD3,U1,ALPH,BETA,CNO,
      T11,T12,T13,T21,T22,T23,T31,TI11,TI12,TI13,TI21,TI22,TI23,TI31,
      TI32,TI33;
    double uround,faccon,fnewt,tolst,xph,firstAcceptedStep;
    double thet,dynold,h,hhfac,theta,dyno,safe,cfac,facr,facl,thqold;
    double hacc,erracc,hold,xold,xend,posneg,hmaxn,hmax,hopt,quot1,quot2;
    bool sameTestValue,startn,first,caljac,reject,gustafssonTest;
    bool last,onlyOneStep;
    fortranVectorF<n> atol,rtol,CONT1,scal;

    fortranVectorF<n> save1,save2,save3;//CONT(N+...) in Hairer.
    MatrixReal Jac;
    int naccpt;
 
    // counters:
    int nit,newt,nmax,nrejct,njac,nstep,ndec,nfonccalled;
  protected:
    Fonct F;
    fortranVectorF<n> Z1,Z2,Z3,F1,F2,F3,Z1T,Z2T,Z3T;
    //! compute (numerically) jacobian of Fonct
    //! \param t current time (IN)
    //! \param  y current value of unknowns (IN)
    //! \param Fy computed values of F(t,y,..) (IN)
    inline void Jacobian(double t,fortranVector y,const fortranVector Fy)
    {
      fortranVectorF<n> R;// move this in class declaration => slower!.
      for(int j=1;j<=n;j++)
	{
	  double ysafe=y(j);
	  double delt=sqrt(uround*MAX(1.e-5,ABS(ysafe))),udelt=1.0e+0/delt;
	  y(j)+=delt;
	  F(t,&y,&R);
	  int i1=Ibegin(j),i2=Iend(j);
#include "Ivdep.hpp"
	  for(int i=i1;i<=i2;i++)
	    Jac(i,j)=(R(i)-Fy(i))*udelt;
	  y(j)=ysafe;
	}
      nfonccalled+=n;
    }
    //! make Newton matrices (factorize them).
    void doNewtonMatrices()
    {
      double fac1=U1/h;
      // compute E1 and E2 matrices, and their decompositions:
      decomr(fac1,Jac);
      double alphn=ALPH/h, betan=BETA/h;
      decomc(alphn,betan,Jac);
      ++ndec;
    }
    //! Linear system solution in Newton iterations.
    void slvrad()
    {
      
      double fac1=U1/h;
      double alphn=ALPH/h, betan=BETA/h;
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)//ZZ coded like in Hairer. optimal??
      	{
      	  double s2=-F2(i),s3=-F3(i);
      	  Z1(i)-=fac1*F1(i);
      	  Z2(i)+=s2*alphn-s3*betan;
      	  CONT1(i)=Z3(i)+s3*alphn+s2*betan;
      	}

      solvereal(Z1,Jac);
      solvecomplex(Z2,CONT1,Jac);
      Z3=CONT1;
    }
    //! Test simplified Newton iterations convergence.
    //! \note return a pair made ok:
    //!  (true or false) depending if the convergence.
    //!  a suggested new value for the time step h.
    std::pair<bool,double> NewtonIterationTest()
    {
      dyno=0.0;
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)
	{
	  double udenom=1./scal(i);
	  dyno+= pow(Z1(i)*udenom,2)+pow(Z2(i)*udenom,2)+pow(Z3(i)*udenom,2);
	}
      dyno=sqrt(dyno/(3*n));
      bool ok=true; double newh=h;
      if(newt>1 && newt<nit)
	{
	  double thq=dyno/dynold;
	  if(newt==2)
	    theta=thq;
	  else
	    theta=sqrt(thq*thqold);
	  thqold=thq;
	  if(theta<0.99)
	    {
	      ok=true;
	      faccon=theta/(1.0-theta);
	      double dyth=faccon*dyno*pow(theta,nit-1-newt)/fnewt;
	      if(dyth>1.0)
		{
		  double qnewt=MAX(1.e-4,MIN(20,dyth));
		  hhfac=.8*pow(qnewt,(-1.0/(4.0+nit-1-newt)));
		  newh*=hhfac;
		  reject=true;
		  last=false;
		  ok=false;
		}
	    }
	  else
	    {
	      // 78 in Hairer.
	      newh=h/2.; hhfac=0.5; reject=true; last=false; ok=false;
	    }
	}
      return make_pair(ok,newh);
    }

    
    //! perform simplified Newton iterations.
    //! \param x  time
    //! \param Y   current value
    //! \note we return the reason for which we break iterations.
    int Newton(double x,fortranVector& Y)
    {
      // simplified Newton iterations:
      newt=0; bool leave=false; int whyDoWeLeave=0;
      faccon=pow(MAX(faccon,uround),0.8);
      theta=ABS(thet);
      do
	{
	
	  CONT1.setsum2(Y,Z1);
	  F(x+C1*h,&CONT1,&Z1);
	  CONT1.setsum2(Y,Z2);
	  F(x+C2*h,&CONT1,&Z2);
	  CONT1.setsum2(Y,Z3); 
	  F(xph,&CONT1,&Z3);
#include "Ivdep.hpp"
	  for(int i=1;i<=n;i++)
	    {
	      double A1=Z1(i);
	      double A2=Z2(i);
	      double A3=Z3(i);
	      Z1T(i)=TI11*A1+TI12*A2+TI13*A3;
	      Z2T(i)=TI21*A1+TI22*A2+TI23*A3;
	      Z3T(i)=TI31*A1+TI32*A2+TI33*A3;
	    }
	  Z1.switchadr(Z1T);Z2.switchadr(Z2T);Z3.switchadr(Z3T);
	  slvrad();
	  nfonccalled+=3;
	  ++newt;
	  //test convergence of Newton iterations:
	  std::pair<bool,double> cv=NewtonIterationTest();
	  if(!cv.first)
	    {
	      seth(x,cv.second);
	      leave=true; whyDoWeLeave=willNotConverge;
	    }
	  else
	    {
	      dynold=MAX(dyno,uround);
	      F1+=Z1; F2+=Z2; F3+=Z3;
	      Z1.lc3(T11,F1,T12,F2,T13,F3);
	      Z2.lc3(T21,F1,T22,F2,T23,F3);
	      Z3.a_x_plus_y(T31,F1,F2);
	      
	      if(faccon*dyno<fnewt)
		{
		  leave=true; whyDoWeLeave=converged;
		}
	    }
	  if(newt==nit)
	    {
	      leave=true; whyDoWeLeave=didNotConverge;
	    }
	}
      while(!leave);
      return whyDoWeLeave;
    }
    //! (re) compute scal
    //! \param Y 
    void doScal(const fortranVector& Y)
    {
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)
	scal(i)=atol(i)+rtol(i)*ABS(Y(i));
    }
    //! error estimation.
    //! \param y0 
    //! \param y
    //! \param t
    //! \param first
    //! \param reject
    double estrad(fortranVectorF<n>& y0,fortranVector& y,
		double t,bool first,bool reject)
    {
      double h1=1.e0/h,n1=1.0e0/n;
      double HEE1=DD1*h1, HEE2=DD2*h1, HEE3=DD3*h1;
      F2.lc3(HEE1,Z1,HEE2,Z2,HEE3,Z3);
      CONT1.sum(F2,y0);
      solvereal(CONT1,Jac);
      //
      double err=0.0;
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)
	err+=pow(CONT1(i)/scal(i),2);
      err=MAX(sqrt(err*n1),1.e-10);
 
      if(err>=1.0)
       if(first||reject)
	 {
	   CONT1.sum(CONT1,y);
	   F(t,&CONT1,&F1);
	   CONT1.sum(F1,F2);
	   solvereal(CONT1,Jac);
	   err=0.0;
	   for(int i=1;i<=n;i++)
	     err+=pow(CONT1(i)/scal(i),2);
	  err=MAX(sqrt(err*n1),1.e-10);
	 }

      return err;
  }
    //! set value of h.
    //! \param x current value of time.
    //! \param newh newvalue of h.
    //! \note can be used to keep an history of h modifications.
    inline void seth(double x,double newh)
    {
      h=newh;
#ifdef LOGRADAU5
      Logger.put(x,nstep,changedH,h);
#endif
    }
  public:
    //! set the error test policy
    //! \param  _sameTestValue == true iff atol and rtol are scalar
    //! \param _atol
    //! \param _rtol
    void setTestPolicy(bool _sameTestValue,double _atol[],double _rtol[])
    {
      sameTestValue=_sameTestValue;
      double expm=2.0/3.0;
      if(sameTestValue)
	{
	  double quot=_atol[0]/_rtol[0];
	  double a=0.1*pow(_rtol[0],expm);
#include "Ivdep.hpp"
	  for(int i=1;i<=n;i++)
	    atol(i)=a*quot;
#include "Ivdep.hpp"
	  for(int i=1;i<=n;i++)
	    rtol(i)=a;
	}
      else
	{
#include "Ivdep.hpp"
	  for(int i=1;i<=n;i++)
	    {
	      double quot=_atol[i-1]/_rtol[i-1];
	      rtol(i)=0.1*pow(_rtol[i-1],expm);
	      atol(i)=rtol(i)*quot;
	    }
	}
      tolst=rtol(1);
      fnewt=MAX(10*uround/tolst,MIN(0.03,sqrt(tolst)));
    }
    //! return a reference to the Jacobian.
    inline MatrixReal& Jacobian(){return Jac;}
    //! return a reference to the right-hand side object function.
    inline Fonct& rhs() {return F;}
#ifdef LOGRADAU5
    //! return a reference to the logger.
    inline logger& Log(){return Logger;}
#endif
    //! constructor
    //! \param _sameTestValue if true error tolerance are scalar
    //! \param _atol absolute tolerance
    //! \param _rtol relative tolerance
    //! \note even if the tolerance are scalar, _atol and _rtol
    //! are arrays.
    Radau5cc(bool _sameTestValue,double _atol[],double _rtol[]): 
	      Matrices<MatrixFull,Hessenberg,n,ninf,nsup>()
    {
      //setTestPolicy(_sameTestValue,_atol,_rtol);
      
      double   sq3=sqrt(3.0);
      SQ6=sqrt(6.0e0);
      C1=(4.-SQ6)/10.;    C2=(4.+SQ6)/10.;
      C1M1=C1-1.;    C2M1=C2-1.;    C1MC2=C1-C2;
      DD1=-(13.+7.*SQ6)/3.;    DD2=(-13.+7.*SQ6)/3.;    DD3=-1./3.;
      U1=(6.e+0+pow(81.e+0,(1.e+0/3.e+0))-pow(9.e+0,(1.e+0/3.e+0)))/30.e+0;
      ALPH=(12.e+0-pow(81.e+0,(1.e+0/3.e+0))+pow(9.e+0,(1.e+0/3.e+0)))/60.e+0;
      BETA=(pow(81.e+0,(1.e+0/3.e+0))+pow(9.e+0,(1.e+0/3.e+0)))*sq3/60.e+0;
      CNO=ALPH*ALPH+BETA*BETA;
      U1=1.0/U1;
      ALPH=ALPH/CNO;    BETA=BETA/CNO;
      T11=9.1232394870892942792e-02;    T12=-0.14125529502095420843;
      T13=-3.0029194105147424492e-02;   T21=0.24171793270710701896;
      T22=0.20412935229379993199;       T23=0.38294211275726193779;
      T31=0.96604818261509293619;       TI11=4.3255798900631553510;
      TI12=0.33919925181580986954;      TI13=0.54177053993587487119;
      TI21=-4.1787185915519047273;      TI22=-0.32768282076106238708;
      TI23=0.47662355450055045196;      TI31=-0.50287263494578687595;
      TI32=2.5719269498556054292;       TI33=-0.59603920482822492497;

      uround=1.e-16;
    
      caljac=true;startn=false; first=true;
      faccon=1.;
      thet=0.001e+0;
      nit=7;
      safe=0.9e+0;
      cfac=safe*(1+2*nit);
    
      facr=1.e+0/8.e+0; facl=5.e+0;
    
      gustafssonTest=true;
    
      nmax=10000;

      setTestPolicy(_sameTestValue,_atol,_rtol);

      //tolst=rtol(1);
      // fnewt=MAX(10*uround/tolst,MIN(0.03,sqrt(tolst)));
      quot1=1.0e+0; quot2=1.2e+0;

      onlyOneStep=false;
    }
    //! destructor
    ~Radau5cc()
    {
    }
    //! operator will perform only one step, and return,whatever the value of
    //! time x.
    //! \param val set onlyOneStep to val
    //! \note default is false, ie operator() will integrate up to _xend
    //! (see operator()).
    inline void performOnlyOneStep(bool val){onlyOneStep=val;}
    inline void setRecomputeJacobianTreshold(double _thet){thet=_thet;}
    inline void setMaxIterationsNewton(int _nit){nit=_nit;}
    //! change safety factor for step size reduction
    //! \param _safe new value
    inline void setSafe(double _safe) {safe=_safe;}
    //! change parameters for step size reduction.
    //! \param _facr
    //! \param _facl
    //! \note new step h is chose so that _facl<= h/hold<=facr
    inline void setFacrFacl(double _facr,double _facl){facr=_facr;facl=_facl;}
    //! perform Gustafsson test? (default yes).
    //! \param val boolean.
    inline void setGustafssonTest(bool val){gustafssonTest=val;}
    //! change the maximum number of allowed steps.
    //! \param  _nmax new value.
    inline void setNmax(int _nmax){nmax=_nmax;}
    //! change stopping criterion for Newton iterations.
    //! \param  _fnewt new value.  
    inline void setFnewt(double _fnewt)
    {
      if(_fnewt<=uround/tolst)
	throw OdesException("Radau5cc:setFnewt, _fnewt<=uround/tolst",
			       "_fnewt,uround,tolst=",_fnewt,uround,tolst);
      fnewt=_fnewt;
    }
    //! get number of steps performed.
    inline int getNstep() const{return nstep;}
    //! get number of steps accepted.
    inline int getNaccpt() const {return naccpt;}
    //! get number of steps rejected.
    inline int getNrejct() const {return nrejct;}
    //! get number of jacobian computed.
    inline int getNJac() const {return njac;}
    //! get number of matrix decompositions.
    inline int getNdec() const {return ndec;}
    //! get last time step used.
    inline double getLastTimeStep() const {return h;}
    //! get number of Newton iterations performed.
    inline int getNewt() const {return newt;}
    //! get number of function call.
    //! \note Z! if jacobian is computed numerically, we have n call for it.
    //! else, an estimation of the cost must take account of the cost of
    //! the computation of the exact Jacobian.
    inline int getNfoncCalled()const {return nfonccalled;}
    //
    inline double getfirstAcceptedStep() const {return firstAcceptedStep;}
    //! make one time step.
    //! \param _h time step
    //! \param x actual time
    //! \param _xend we integrate from x to _xend 
    //! \param y[] current unknowns values.
    void operator()(double _h,double& x,double _xend,double y[])
    {
#ifdef LOGRADAU5
      Logger.clear();
#endif
      //initialize counters:
      naccpt=0; nrejct=0; njac=0;nfonccalled=0;
      //-------------------
      bool redonewton=true; 
      caljac=true;startn=false; first=true;
      reject=false,last=false;
      cfac=safe*(1+2*nit); facr=1.e+0/8.e+0; facl=5.e+0;
      int reason=0;
      nstep=0; ndec=0;
      //h=_h;
      seth(x,_h);
      fortranVectorF<n> y0;
      fortranVector Y(y,n);
      xend=_xend;
      posneg=SIGN(xend-x);
      hmax=xend-x;
      hmaxn=MIN(ABS(hmax),ABS(xend-x));
      if(ABS(h)<10.*uround) //h=1.e-6;
	seth(x,1.e-6);
      seth(x,MIN(ABS(h),hmaxn)*posneg);
      hold=h; xold=x;
      if((x+h*1.0001-xend)*posneg>0.0)
	{
	  seth(x,xend-x);
	  last=true;
	}
      else
      hopt=h; xold=x;

      F(x,y,&y0);
      ++nfonccalled;
      hhfac=h;

      // compute Jacobian:
      if(ComputeJacobianNumerically)
	Jacobian(x,Y,y0);
      else
	F.Jacobian(x,Y,y0,Jac);
      ++njac; caljac=true; calhes=true;

      doNewtonMatrices();
      doScal(Y);

      for(;;)
	{
	  if(nstep>nmax)
	    throw OdesException("Radau5cc() nstep=",nstep,"too large, nmax=",
				   nmax);

	  do{//20 in Hairer.
	    ++nstep;
	    // starting values for newton iteration:
	    xph=x+h;
	    if(startn||first)
	      {
		Z1=0.0; Z2=0.0; Z3=0.0;
		F1=0.0; F2=0.0; F3=0.0;
	      }
	    else
	      {
		double c3q=h/hold,c1q=C1*c3q,c2q=C2*c3q;
		double c1q_c2m1=c1q-C2M1,c1q_c1m1=c1q-C1M1;
		double c2q_c2m1=c2q-C2M1,c2q_c1m1=c2q-C1M1;
		double c3q_c2m1=c3q-C2M1,c3q_c1m1=c3q-C1M1;
#include "Ivdep.hpp"
		for(int i=1;i<=n;i++)
		  {
		    Z1(i)=
		      c1q*(save1(i)+c1q_c2m1*(save2(i)+c1q_c1m1*save3(i)));
		    Z2(i)=
		      c2q*(save1(i)+c2q_c2m1*(save2(i)+c2q_c1m1*save3(i)));
		    Z3(i)=
		      c3q*(save1(i)+c3q_c2m1*(save2(i)+c3q_c1m1*save3(i)));
		  }


		F1.lc3(TI11,Z1,TI12,Z2,TI13,Z3);
		F2.lc3(TI21,Z1,TI22,Z2,TI23,Z3);
		F3.lc3(TI31,Z1,TI32,Z2,TI33,Z3);

	      }
	    // perform simplified Newton iterations.
	    reason=Newton(x,Y);

	    //test how did we exited the Newton loop.
	    if(reason==converged)
	      redonewton=false;
	    else
	      {
		if(reason==didNotConverge)
		  {
		    //h/=2.; 
		    seth(x,h*0.5);
		    hhfac=0.5;
		    if(!caljac) 
		      {
			if(ComputeJacobianNumerically)
			  Jacobian(x,Y,y0);
			else
			  F.Jacobian(x,Y,y0,Jac);
			++njac; caljac=true; calhes=true;
		      }
		    doNewtonMatrices();
		    doScal(Y);
		    newt=0;
#ifdef LOGRADAU5
		    Logger.put(x,nstep,NewtonFailed,h);
#endif
		  }
		else if(reason==willNotConverge)
		  {
		    newt=0;
		    redonewton=true;
		    if(!caljac) 
		      {
			if(ComputeJacobianNumerically)
			  Jacobian(x,Y,y0);
			else
			  F.Jacobian(x,Y,y0,Jac);
			++njac; caljac=true; calhes=true;
		      }
		    doNewtonMatrices();
		    doScal(Y);
#ifdef LOGRADAU5
		    Logger.put(x,nstep,NewtonWillNotConverge,h);
#endif		  
		  }
	      }
	  }
	  while(redonewton);
	  // error estimation
	  double err=estrad(y0,Y,x,first,reject);
	  // new h.
	  double fac=MIN(safe,cfac/(newt+2*nit));
	  double quot=MAX(facr,MIN(facl,sqrt(sqrt(err))/fac));
	  double hnew=h/quot;
	  if(err<1.0)
	    {
	      if(onlyOneStep) {Y+=Z3;break;}
	      //step accepted
	      if(naccpt==0) firstAcceptedStep=h;
	      first=false;++naccpt;
	      //Gustafsson's step control:
	      if(gustafssonTest)
		{
		  if(naccpt>1)
		    {
		      double facgus= (hacc/h)*sqrt(sqrt((err*err/erracc)))/safe;
		      facgus=MAX(facr,MIN(facl,facgus));
		      quot=MAX(quot,facgus);
		      hnew=h/quot;
		    }
		  hacc=h;
		  erracc=MAX(0.01,err);
		}
	      xold=x,hold=h;
	      x=xph;
	      Y+=Z3;
	      
	      double uC2M1=1./C2M1,uC2=1./C2;
#include "Ivdep.hpp"
	      for(int i=1;i<=n;i++)
		save1(i)=(Z2(i)-Z3(i))*uC2M1;
	      double uC1MC2=1./C1MC2,uC1M1=1./C1M1;
#include "Ivdep.hpp"
	      for(int i=1;i<=n;i++)
		save2(i)=((Z1(i)-Z2(i))*uC1MC2-save1(i))*uC1M1;
	      double aaa=(uC1MC2-uC1M1)*uC2,bbb=uC1MC2*uC2;
#include "Ivdep.hpp"
	      for(int i=1;i<=n;i++)
		save3(i)=save2(i)-aaa*Z1(i)+bbb*Z2(i); 

	      doScal(Y);
	      caljac=false;
	      if(last||0.1*ABS(xend-x)<=ABS(x)*uround)
		{
		  seth(x,hopt);
		  break;
		}
	      F(x,y,&y0);
	      ++nfonccalled;
	      hnew=posneg*MIN(ABS(hnew),hmaxn);
	      hopt=MIN(h,hnew);
	      if(reject)
		hnew=posneg*MIN(ABS(hnew),ABS(h));
	      reject=false;

	      if((x+hnew/quot1-xend)*posneg>=0.0)
		{
		  seth(x,xend-x);
		  last=true;
		  hhfac=h;
		  if(theta>thet)
		    {
		      if(ComputeJacobianNumerically)
			Jacobian(x,Y,y0);
		      else
			F.Jacobian(x,Y,y0,Jac);
		      ++njac; caljac=true; calhes=true;
		    }
		  doNewtonMatrices();

		}
	      else
		{
		  double qt=hnew/h;
		  hhfac=h;
		  if(!(theta<=thet && qt>=quot1 && qt<=quot2))
		    {
		      if(theta>thet)
			{
			  if(ComputeJacobianNumerically)
			    Jacobian(x,Y,y0);
			  else
			    F.Jacobian(x,Y,y0,Jac);
			  ++njac; caljac=true; calhes=true;
			}
		      seth(x,hnew);
		      doNewtonMatrices();
		    }
		}
	    }
	  else
	    {
	      //step rejected
#ifdef LOGRADAU5
	      Logger.put(x,nstep,rejectedStep,h);
#endif
	      reject=true; last=false;
	      if(first)
		{
		  h*=0.1; hhfac=0.1;
#ifdef LOGRADAU5
		  Logger.put(x,nstep,changedH,h);
#endif
		}
	      else
		{
		  hhfac=hnew/h;
		  h=hnew;
#ifdef LOGRADAU5
		  Logger.put(x,nstep,changedH,h);
#endif
		}
	      if(naccpt) ++nrejct;
	      if(!caljac)
		{
		  if(ComputeJacobianNumerically)
		    Jacobian(x,Y,y0);
		  else
		    F.Jacobian(x,Y,y0,Jac);
		  ++njac; caljac=true; calhes=true;
		}
	      doNewtonMatrices();
	    }
	  
	}//end for loop.
    }
  };
}
#endif
