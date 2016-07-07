#ifndef RODASCC__H
#define RODASCC__H
#include "OdesException.hpp"
#include "fortranArray.hpp"
#include "fortranVector.hpp"
#include "protos_lapack.hpp"
#include "Matrices.hpp"
#include "compat.hpp"
#include <cmath>
#include <utility>
#ifdef LOGRODAS
#include "logger.hpp"
#endif
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
#define SIGN(a) (((a) >= (0.0)) ? (1.) : (-1.))
namespace odes
{
  /////////////////////////////////////////////////////////////////////
  /// Class Rodascc
  ///
  /// A C++ implementation of (most part of) Hairer & Wanner Rodas 
  /// program.
  /// \brief A C++ implementation of Rodas.
  ////////////////////////////////////////////////////////////////////
  template<class Fonct> class Rodascc: 
    private Matrices<(Fonct::n-Fonct::nsub)==1&&(Fonct::n-Fonct::nsup)==1,
		     Fonct::Hessenberg,Fonct::n,Fonct::nsub,Fonct::nsup>
  {
    static const int n=Fonct::n;
    static const int ninf=Fonct::nsub;
    static const int nsup=Fonct::nsup;
    static const bool MatrixFull=(n-ninf)==1&&(n-nsup)==1;
    static const bool Hessenberg=Fonct::Hessenberg;
    static const bool ComputeJacobianNumerically=
      Fonct::ComputeJacobianNumerically;
    static const int meth=Fonct::method;// which method do we use?
    static const bool autonms=Fonct::autonomous;// autonomous problem?
    static const bool use_DF_t=Fonct::use_DF_t;// analytical or not.

    typedef typename Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::MatrixReal
    MatrixReal;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::Ibegin;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::Iend;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::calhes;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::decomr;
    using Matrices<MatrixFull,Hessenberg,n,ninf,nsup>::solvereal;
#ifdef LOGRODAS
    logger Logger;
#endif
    // We try to keep most of the notation of Hairer & Wanner.
    double A21,A31,A32,A41,A42,A43,A51,A52,A53,A54,
      C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,
      C62,C63,C64,C65,GAMMA,C2,C3,C4,D1,D2,D3,D4,
      D21,D22,D23,D24,D25,D31,D32,D33,D34,D35,HD1,HD2,HD3,HD4;
    double uround,faccon,fac,fac1,fac2,h,safe;
    double xend,posneg,hmax,xold,hacc,erracc,hmaxn,hopt,firstAcceptedStep;

    bool sameTestValue,caljac,startn,first,gustafssonTest,reject;
    bool last,onlyOneStep;

    fortranVectorF<n> AK,AK1,AK2,AK3,AK4,AK5,AK6,DY,DY1,YNEW,FX,
      atol,rtol,CONT;
 
    MatrixReal Jac;
    // counters:
    int naccpt,nfcn,nsol,newt,nmax,nrejct,njac,nstep,ndec;
 
    //define coefficients:
    void rocoe()
    {
      switch(meth)
	{
	case 1:
	  C2=0.386e0;C3=0.21e0 ;C4=0.63e0;
	  D1= 0.2500000000000000e+00;  D2=-0.1043000000000000e+00;
	  D3= 0.1035000000000000e+00;  D4=-0.3620000000000023e-01;
	  A21= 0.1544000000000000e+01; A31= 0.9466785280815826e+00;
	  A32= 0.2557011698983284e+00; A41= 0.3314825187068521e+01;
	  A42= 0.2896124015972201e+01; A43= 0.9986419139977817e+00;
	  A51= 0.1221224509226641e+01; A52= 0.6019134481288629e+01;
	  A53= 0.1253708332932087e+02; A54=-0.6878860361058950e+00;
	  C21=-0.5668800000000000e+01; C31=-0.2430093356833875e+01;
	  C32=-0.2063599157091915e+00; C41=-0.1073529058151375e+00;
	  C42=-0.9594562251023355e+01; C43=-0.2047028614809616e+02;
	  C51= 0.7496443313967647e+01; C52=-0.1024680431464352e+02;
	  C53=-0.3399990352819905e+02; C54= 0.1170890893206160e+02;
	  C61= 0.8083246795921522e+01; C62=-0.7981132988064893e+01;
	  C63=-0.3152159432874371e+02; C64= 0.1631930543123136e+02;
	  C65=-0.6058818238834054e+01;
	  GAMMA= 0.2500000000000000e+00;
	  D21= 0.1012623508344586e+02; D22=-0.7487995877610167e+01;
	  D23=-0.3480091861555747e+02; D24=-0.7992771707568823e+01;
	  D25= 0.1025137723295662e+01; D31=-0.6762803392801253e+00;
	  D32= 0.6087714651680015e+01; D33= 0.1643084320892478e+02;
	  D34= 0.2476722511418386e+02; D35=-0.6594389125716872e+01;
	  break;
	case 2:
	  C2=0.3507221;	  C3=0.2557041;	  C4=0.6817790;

	  D1= 0.2500000000000000e+00;  D2=-0.6902209999999998e-01;
	  D3=-0.9671999999999459e-03;  D4=-0.8797900000000025e-01;
	  A21= 0.1402888400000000e+01; A31= 0.6581212688557198e+00;
	  A32=-0.1320936088384301e+01; A41= 0.7131197445744498e+01;
	  A42= 0.1602964143958207e+02; A43=-0.5561572550509766e+01;
	  A51= 0.2273885722420363e+02; A52= 0.6738147284535289e+02;
	  A53=-0.3121877493038560e+02; A54= 0.7285641833203814e+00;
	  C21=-0.5104353600000000e+01; C31=-0.2899967805418783e+01;
	  C32= 0.4040399359702244e+01; C41=-0.3264449927841361e+02;
	  C42=-0.9935311008728094e+02; C43= 0.4999119122405989e+02;
	  C51=-0.7646023087151691e+02; C52=-0.2785942120829058e+03;
	  C53= 0.1539294840910643e+03; C54= 0.1097101866258358e+02;
	  C61=-0.7629701586804983e+02; C62=-0.2942795630511232e+03;
	  C63= 0.1620029695867566e+03; C64= 0.2365166903095270e+02;
	  C65=-0.7652977706771382e+01; GAMMA= 0.2500000000000000e+00  ;
	  
	  D21=-0.3871940424117216e+02; D22=-0.1358025833007622e+03;
	  D23= 0.6451068857505875e+02; D24=-0.4192663174613162e+01;
	  D25=-0.2531932050335060e+01; D31=-0.1499268484949843e+02;
	  D32=-0.7630242396627033e+02; D33= 0.5865928432851416e+02;
	  D34= 0.1661359034616402e+02; D35=-0.6758691794084156e+00;
	  break;
	case 3:
	  GAMMA = 0.25;	C2=3.0*GAMMA; C3=0.21; C4=0.63;

	  D1= 0.2500000000000000e+00; D2=-0.5000000000000000e+00;
	  D3=-0.2350400000000000e-01; D4=-0.3620000000000000e-01;
	  A21= 0.3000000000000000e+01;A31= 0.1831036793486759e+01;
	  A32= 0.4955183967433795e+00;A41= 0.2304376582692669e+01;
	  A42=-0.5249275245743001e-01;A43=-0.1176798761832782e+01;
	  A51=-0.7170454962423024e+01;A52=-0.4741636671481785e+01;
	  A53=-0.1631002631330971e+02;A54=-0.1062004044111401e+01;
	  C21=-0.1200000000000000e+02;C31=-0.8791795173947035e+01;
	  C32=-0.2207865586973518e+01;C41= 0.1081793056857153e+02;
	  C42= 0.6780270611428266e+01;C43= 0.1953485944642410e+02;
	  C51= 0.3419095006749676e+02;C52= 0.1549671153725963e+02;
	  C53= 0.5474760875964130e+02;C54= 0.1416005392148534e+02;
	  C61= 0.3462605830930532e+02;C62= 0.1530084976114473e+02;
	  C63= 0.5699955578662667e+02;C64= 0.1840807009793095e+02;
	  C65=-0.5714285714285717e+01;

	  D21= 0.2509876703708589e+02;D22= 0.1162013104361867e+02;
	  D23= 0.2849148307714626e+02;D24=-0.5664021568594133e+01;
	  D25= 0.0000000000000000e+00;D31= 0.1638054557396973e+01;
	  D32=-0.7373619806678748e+00;D33= 0.8477918219238990e+01;
	  D34= 0.1599253148779520e+02;D35=-0.1882352941176471e+01;
	  break;
	default:
	    throw OdesException("Rodascc, meth=",meth);
	  
	}

    }
  protected:
    Fonct F;
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
	  int id=Ibegin(j),iend=Iend(j);
#include "Ivdep.hpp"
	  for(int i=id;i<=iend;i++)
	    Jac(i,j)=(R(i)-Fy(i))*udelt;
	  y(j)=ysafe;
	}
    }
    //compute derivative of F/x
    //! \param x  we solve h'=F(x,y)
    //! \param Y the result (OUT).
    inline void DF_X(double x,const fortranVectorF<n>& Y)
    {
      if(use_DF_t)
	F.DF_t(x,Y,FX);
      else
	{
	  double delt=sqrt(uround*MAX(1.e-5,ABS(x)));
	  F(x+delt,&Y,&AK1);
	  double delt1=1./delt;
	  FX.lc2(delt1,AK1,-delt1,DY1);
	}
    }
    //! Linear system solution in Newton iterations.
    void slvrod(double hd,fortranVectorF<n>& A,
		fortranVectorF<n>& DY,fortranVectorF<n>& FX,
		fortranVectorF<n>& Y,bool addYnew)
    {
      if(autonms)
	A=DY;
      else
	A.a_x_plus_y(hd,FX,DY);
      if(addYnew) A+=Y;
      solvereal(A,Jac);
    }
    //! set value of h.
    //! \param x current value of time.
    //! \param newh newvalue of h.
    //! \note can be used to keep an history of h modifications.
    inline void seth(double x,double newh)
    {
      h=newh;
#ifdef LOGRODAS
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
      if(sameTestValue)
#include "Ivdep.hpp"
	for(int i=1;i<=n;i++)
	  {
	    atol(i)=_atol[0];
	    rtol(i)=_rtol[0];
	  }
      else
	{
#include "Ivdep.hpp"
	  for(int i=1;i<=n;i++)
	    {
	      rtol(i)=_rtol[i];
	      atol(i)=_atol[i];
	    }
	}
    }
    //! return a reference to the Jacobian.
    inline MatrixReal& Jacobian(){return Jac;}
    //! return a reference to the right-hand side object function.
    inline Fonct& rhs() {return F;}
#ifdef LOGRODAS
    //! return a reference to the logger.
    inline logger& Log(){return Logger;}
#endif
    //! constructor
    //! \param _sameTestValue if true error tolerance are scalar
    //! \param _atol absolute tolerance (IN)
    //! \param _rtol relative tolerance (IN)
    //! \note even if the tolerance are scalar, _atol and _rtol
    //! are arrays.
    Rodascc(bool _sameTestValue,double _atol[],double _rtol[]): 
	      Matrices<MatrixFull,Hessenberg,n,ninf,nsup>()
    {
      
      uround=1.e-16;
    
      caljac=true;startn=false; first=true;
      faccon=1.;

      safe=0.9e+0;
      fac1=5.; fac2=1./6.;
      gustafssonTest=true;
    
      nmax=10000;

      setTestPolicy(_sameTestValue,_atol,_rtol);

      rocoe();
      
      onlyOneStep=false;
    }
    //! destructor
    ~Rodascc()
    {
    }
    //! operator will perform only one step, and return,whatever the value of
    //! time x.
    //! \param val set onlyOneStep to val
    //! \note default is false, ie operator() will integrate up to _xend
    //! (see operator()).
    inline void performOnlyOneStep(bool val){onlyOneStep=val;}
    //! change safety factor for step size reduction
    //! \param _safe new value
    inline void setSafe(double _safe) {safe=_safe;}
    //! perform Gustafsson test? (default yes).
    //! \param val boolean.
    inline void setGustafssonTest(bool val){gustafssonTest=val;}
    //! change the maximum number of allowed steps.
    //! \param  _nmax new value.
    inline void setNmax(int _nmax){nmax=_nmax;}
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
    //
    inline double getfirstAcceptedStep() const {return firstAcceptedStep;}
    //! get last time step used.
    inline double getLastTimeStep() const {return h;}
    //! make one time step.
    //! \param _h time step
    //! \param x actual time
    //! \param _xend we integrate from x to _xend 
    //! \param y[] current unknowns values.
    void operator()(double _h,double& x,double _xend,double y[])
    {
#ifdef LOGRODAS
      Logger.clear();
#endif
      double HD1,HD2,HD3,HD4;
      //initialize counters:
      naccpt=0; nrejct=0; njac=0;
      //-------------------
      startn=false; first=true;
      reject=false,last=false;
      nstep=0; ndec=0;
      nsol=0; nfcn=0;

      seth(x,_h);
      fortranVectorF<n> y0;
      fortranVectorF<n> Y;
      Y=y;
      xend=_xend;
      posneg=SIGN(xend-x);
      hmax=xend-x;
      hmaxn=MIN(ABS(hmax),ABS(xend-x));
      while(!last)
	{
	  if (0.1e0*ABS(h)<=ABS(x)*uround)
	    throw OdesException("Rodascc():, h too small",h);
	  if(ABS(h)<10.*uround) //h=1.e-6;
	    seth(x,1.e-6);
	  seth(x,MIN(ABS(h),hmaxn)*posneg);
	  xold=x;
	  if((x+h*1.0001-xend)*posneg>0.0)
	    {
	      seth(x,xend-x);
	      last=true;
	    }
	  else
	    hopt=h; xold=x;

	  F(x,&Y,&DY1);
	  
	  //compute Jacobian:
	  if(ComputeJacobianNumerically)
	    Jacobian(x,Y,DY1);
	  else
	    F.Jacobian(x,Y,DY1,Jac);
	  ++njac; caljac=true; calhes=true;
	  
	  // in the non autonomous case, compute derivative /x.

	  if(!autonms) DF_X(x,Y);

	  do
	    {
	      fac=1./(h*GAMMA);
	     
	      decomr(fac,Jac);
	      ++ndec;
	      //prepare computation of the 5 stages:
	      double HC21=C21/h,HC31=C31/h, HC32=C32/h, HC41=C41/h, HC42=C42/h,
		HC43=C43/h, HC51=C51/h, HC52=C52/h, HC53=C53/h, HC54=C54/h, 
		HC61=C61/h, HC62=C62/h, HC63=C63/h, HC64=C64/h, HC65=C65/h;
	      if(autonms)
		{
		  HD1=0.0; HD2=0.0; HD3=0.0; HD4=0.0;
		}
	      else
		{
		  HD1=h*D1; HD2=h*D2; HD3=h*D3; HD4=h*D4;
		}

	      //the stages:
	      // 1:
	      slvrod(HD1,AK1,DY1,FX,YNEW,false);
	      YNEW.a_x_plus_y(A21,AK1,Y);
	      //
	      //2:
	      F(x+C2*h,&YNEW,&DY);
	      YNEW.a_Y(HC21,AK1);
	      slvrod(HD2,AK2,DY,FX,YNEW,true);
	      YNEW.x_lc2(Y,A31,AK1,A32,AK2);
	      //
	      // 3:
	      F(x+C3*h,&YNEW,&DY);
	      YNEW.lc2(HC31,AK1,HC32,AK2);
	      slvrod(HD3,AK3,DY,FX,YNEW,true);
	      YNEW.x_lc3(Y,A41,AK1,A42,AK2,A43,AK3);
	     
	      //
	      // 4:
	      F(x+C4*h,&YNEW,&DY);
	      YNEW.lc3(HC41,AK1,HC42,AK2,HC43,AK3);
	      slvrod(HD4,AK4,DY,FX,YNEW,true);
	      YNEW.x_lc4(Y,A51,AK1,A52,AK2,A53,AK3,A54,AK4);
      
	      //5:
	      F(x+h,&YNEW,&DY);
	      AK6.lc4(HC52,AK2,HC54,AK4,HC51,AK1,HC53,AK3);
	      slvrod(0.0,AK5,DY,FX,AK6,true);

	      // embedded solution:
	      YNEW+=AK5;
	      F(x+h,&YNEW,&DY);
	      CONT.lc5(HC61,AK1,HC62,AK2,HC65,AK5,HC64,AK4,HC63,AK3);
	      slvrod(0.0,AK6,DY,FX,CONT,true);
	      
	      // new solution
	      YNEW+=AK6;
	      nsol+=6; nfcn+=5;
	      //error estimation:
	      ++nstep; 
	      double err=0.0;
	      for(int i=1;i<=n;i++)
		{
		  double q=atol(i)+rtol(i)*MAX(ABS(Y(i)),ABS(YNEW(i)));
		  err+=pow(AK6(i)/q,2);
		}
	      err=sqrt(err/n);
	      
	      // computation of hnew. 
	      // We require: .2<=hnew/h<=6.
	      fac=MAX(fac2,MIN(fac1,sqrt(sqrt(err))/safe));
	      double hnew=h/fac;
	      //Is error small enough ?
	      if(err<=1.0)
		{
		  if(naccpt==0) firstAcceptedStep=h;
		  ++naccpt;
		  if(gustafssonTest)
		    {
		      if(naccpt>1)
			{
			  double facgus=
			    (hacc/h)* sqrt(sqrt(err*err/erracc))/safe;
			  facgus=MAX(fac2,MIN(fac1,facgus));
			  hnew=h/MAX(fac,facgus);
			}
		      hacc=h;
		      erracc=MAX(1.e-02,err);
		    }
		  Y=YNEW;
		  xold=x;
		  x+=h;
		  if(ABS(hnew)>hmaxn) hnew=posneg*hmaxn;
		  if(reject) hnew=posneg*MIN(ABS(hnew),ABS(h));
		  reject=false;
		  h=hnew;
		}
	      else
		{
		  reject=true;
		  last=false;
		  h=hnew;
		  if(naccpt>1) ++nrejct;
		}
	    }
	  while(reject);
	  if(nstep>nmax)
	    throw OdesException("Rodas(): nstep too large:",nstep);
	  if(onlyOneStep) break;
	}
      
      Y.put(y);
    }
  };
}
#endif
