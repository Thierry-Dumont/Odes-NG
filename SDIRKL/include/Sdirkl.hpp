#ifndef Sdirkl__h
#define Sdirkl__h
#include "AllocateDestroyVector.hpp"
#include "GenericException.hpp"
#include <iostream>
using namespace std;
namespace odes
{
///////////////////////////////////////////////////////////////////////////
/// General SDIRK Method for Linear problems.
///////////////////////////////////////////////////////////////////////////
/// @tparam OpLin the Linear operator.
/// @tparam RKMethod the Runge-Kutta method (see for example RKMethod1.hpp)
///////////////////////////////////////////////////////////////////////////
template<class RKMethod,class OpLin> class Sdirkl
{
  static const int nstage=RKMethod::nstage;//!< number of stages.
  RKMethod RK;
  OpLin Op;
  const int n; //!< size of the problems (number of unknowns).
  const double dt; //!< time step
  const bool affine; //!< true in the affine case, when solve du/dt=Au and B!=0.
  double *z; //!< to store all the steps.
  double AD[nstage][nstage];//!< coefficients of RK method/ diagonal.
  double dloc[nstage],aloc[nstage];
  double *Zused[nstage],*Zi[nstage],*Zcomb[nstage+1];
  double *Work,*Work1;
  double diag;
  int lused;//!< how many non zero b's in the RK formula.
  //! copy x into y.
  //! !\param x
  //! !\param y
  void copy(double x[],double y[])
  {
#include "Ivdep.hpp"
    for(int i=0;i<n;i++) y[i]=x[i];
  }
  //! Linear combinaison of nv vectors.
  //!\param nv
  //!\param c[] coefficients.
  //!\param v[] adresses of vectors.
  //!\param out result: \f$ out =  \sum_{i=1}^{nv} c_i v_i.\f$
  void cl(int nv, double c[], double* v[], double out[])
  {
    double *v0=v[0];
    for(int i=0;i<n;i++)
      {
	out[i]=c[0]*v0[i];
#include "Ivdep.hpp"
	for(int j=1;j<nv;j++)
	    out[i]+=v[j][i]*c[j];
      }
  }
  //! Add a linear combinaison of nv vectors to a vector.
  //!\param nv
  //!\param c[] coefficients.
  //!\param v[] adresses of arrays to combine.
  //!\param out result: \f$ out +=  \sum_{i=1}^{nv} c_i v_i.\f$
  void pluscl(int nv, double c[], double* v[], double out[])
  {
    for(int i=0;i<n;i++)
#include "Ivdep.hpp"
      for(int j=0;j<nv;j++)
       	out[i]+=c[j]*v[j][i];
  }
  //! out= a*b[]+c[]
  //! \param a
  //! \param b
  //! \param c
  //! \param out
  void plusa(double a,double b[],double c[],double out[])
  {
#include "Ivdep.hpp"
    for(int i=0;i<n;i++)
      out[i]=a*b[i]+c[i];
  }
  //! z-=b.
  //!\param z
  //!\param b
  void minus(double z[], double b[])
  {
#include "Ivdep.hpp"
    for(int i=0;i<n;i++)
      z[i]-=b[i];
  }
  //! z+=b
  //! \param z
  //! \param b
  void plus(double z[], double b[])
  {
#include "Ivdep.hpp"
    for(int i=0;i<n;i++)
      z[i]+=b[i];
  }
  //! keep only non zero coefficients for the final step.
  //! \note return the number of non zero coefficients.
  int doNew()
  {
    int ret=0;
    for(int i=0;i<nstage;i++)
      if(RK.D[i]!=0.0)
      {
	  dloc[ret]=RK.D[i]; Zused[ret]=Zi[i]; ++ret;
      }
 
    return ret;
  }
public:
  //!\ Constructor.
  //!\param _n number of unknowns.
  //!\param _dt time step
  //!\param _affine true iff we solve du/dt=Au+B and B!=0.
  Sdirkl(int _n, double _dt,bool _affine=false): n(_n),dt(_dt),affine(_affine)
  {
    if(affine)
      z=new double[n*(nstage+2)];
    else
      z=new double[n*(nstage+1)];
    for(int i=0;i<nstage;i++)
      {
	Zi[i]=z+i*n;
	Zcomb[i+1]=z+i*n;
      }
    Work=z+n*nstage;
    if(affine) Work1=Work+n;
    diag=RK.A[0][0];
    Op.init(n,dt*diag);
    for(int i=0;i<nstage;i++)
#include "Ivdep.hpp"
      for(int j=0;j<=i;j++)
	AD[i][j]=RK.A[i][j]/diag;
    lused=doNew();
  }
  //! destructor.
  ~Sdirkl()
  {
    delete[] z; //delete[] Work;
  }
  //! perform one step.
  //!\param inout initial values (in); computed values (out).
  void step(double inout[])
  {
    if(affine)
      throw GenericException("Sdirkl step(inout) called in the affine case");
    Zcomb[0]=inout;
    for(int stage=0;stage<nstage;stage++)
      {
	aloc[0]=RK.S[stage]/diag;
#include "Ivdep.hpp"
	for(int j=0;j<stage;j++)
	  {
	    aloc[j+1]=AD[stage][j];
	    Zcomb[j+1]=Zi[j];
	  }
	cl(stage+1,aloc,Zcomb,Work);
	Op.apply(Work,Zi[stage]);
	minus(Zi[stage],Work);
      }

    pluscl(lused,dloc,Zused,inout);
  }
  //! perform one step, with a forcing term.
  //!\param inout initial values (in); computed values (out).
  //!\param phi forcing term (indepedent of t).
  void step(double inout[], double phi[])
  {
   if(!affine)
     throw GenericException("Sdirkl step(inout,phi) called",
			    "in the linear case");
    Zcomb[0]=inout;
    for(int stage=0;stage<nstage;stage++)
      {
	aloc[0]=RK.S[stage]/diag;
#include "Ivdep.hpp"
	for(int j=0;j<stage;j++)
	  {
	    aloc[j+1]=AD[stage][j];
	    Zcomb[j+1]=Zi[j];
	  }
	cl(stage+1,aloc,Zcomb,Work1);
	plusa(RK.S[stage]*dt,phi,Work1,Work);
	Op.apply(Work,Zi[stage]);
	minus(Zi[stage],Work1);
      }

    pluscl(lused,dloc,Zused,inout);
  }
};
}
#endif
