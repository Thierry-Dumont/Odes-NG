#include <iostream>
#include "Sdirkl.hpp"
#include "RKMethod1.hpp"
#include "RKMethod2.hpp"
#include "RKRS.hpp"
#include "EulerImp.hpp"
#include "CrN.hpp"
#include "Lapl1d.hpp"
#include <fstream>
#include <cmath>
#include <string>
using namespace std;
using namespace odes;
class exact{
  const int l=25;
  const int n;
  const double pi;
  const double h;
public:
  exact(int _n):n(_n),pi(4.*atan(1.0)),h(1./((double) (n-1)))
  {
  }
  ~exact(){}
  void cinit(double x[])
  {
    double cpi=pi*l/((double) (n-1));
    for(int i=0;i<n;i++)
      x[i]=cos(cpi*i);
 
  }
  void initsecm(double f[])
  {
    cinit(f);
  }
  void exactsolAtt(double t,double e[],double beta)
  {
    cinit(e);
    double elambda=exp(eigenvalue()*t);
    double lambda=eigenvalue();
    double c=elambda*(1+beta/lambda)-beta/lambda;
    for(int i=0;i<n;i++)
      e[i]*=c;
  }
  double eigenvalue() const {return -2.*(1-cos(pi*l/((double)(n-1))))/(h*h);}
};
double normdiff(int n,double x[],double y[])
{
  double ret=0.0;
  for(int i=0;i<n;i++)
    ret+=abs(x[i]-y[i]);
  return ret;
}
int main()
{
  // RK method:
  //typedef RKMethod1 RK;
  //typedef RKMethod2 RK;
  //typedef RKRS RK;
  typedef CrN RK;
  //typedef EulerImp RK;
  

  // Problem:
  typedef Lapl1d Prob;

  //
  int n=1000;//size
 
  double *x=new double[n];
  double *sol=new double[n];
  double *force=new double[n];
  ofstream f;

  double Tinteg=0.01;
  double beta=1.0;

  RK R; //just to get the name!
  f.open(R.name.c_str());
  exact ExactSol(n);
  ExactSol.exactsolAtt(Tinteg,sol,beta);
 

  for(int l=0;l<8;l++)
    {
      ExactSol.cinit(x);
      ExactSol.cinit(force);
      for(int i=0;i<n;i++)
	force[i]*=beta;
      int np=pow(2,l);
      double dt=Tinteg/(double) np;
      Sdirkl<RK,Prob> S(n,dt,true);
      for(int s=0;s<np;s++)
	S.step(x,force);
	//S.step(x);
      double diff=normdiff(n,x,sol);
      cout<<endl<<"step: "<<dt<<" "<<diff<<endl;
      f<<log(dt)<<" "<<log(diff)<<endl;
 
    }
  f.close();
  cout<<"end"<<endl;
}
