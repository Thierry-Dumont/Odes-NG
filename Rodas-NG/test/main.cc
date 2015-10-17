//#define LOGRADAU5
#include <iostream>
#include "AllocateDestroyVector.hpp"
#include "Rodascc.hpp"
#include "AvcF.hpp"
#include "Oregonator.hpp"
#include "KPP.hpp"
#include "E5.hpp"
using namespace std;
int main()
{
  //---- define here the problem you want to treat:
  typedef Oregonator Fonc;
  //typedef KPP Fonc;
  //typedef AvcF Fonc;
  //typedef E5 Fonc;
  //----------------------------------------------
  cout.precision(17);

  static const int n=Fonc::n;
  double atol[n],rtol[n]; atol[0]=1.e-5; rtol[0]=1.e-5;

  Rodascc<Fonc> Rod(true,atol,rtol);
  double h,t,xend;
  int nloops;
  double *y=allocDoubleArray(n);
  //cout<<"initial time step?"; cin>>h;
  //cout<<"integration time?";  cin>>xend;
  h=1.; xend=100.;nloops=1000;
  
  for(int i=0;i<nloops;i++)
    {
      Rod.rhs().init(y);
      t=0.0;
      Rod(h,t,xend,&(y[0]));
    }

  cout<<"nstep= "<<Rod.getNstep()<<" njac= "<<Rod.getNJac()<<
    " naccpt= "<<Rod.getNaccpt()<<" nreject= "<<Rod.getNrejct()<<
    " ndec= "<<Rod.getNdec()<<endl;
  cout<<"first accepted step: "<<Rod.getfirstAcceptedStep()<<endl;
  cout<<"y= ";
  for(int i=0;i<n;i++)
    cout<<y[i]<<" "; cout<<endl;
  cout<<"h= "<<h<<"  t= "<<t<<endl;
#ifdef LOGRADAU5
  Rod.Log().print();
#endif
  cout<<"ok"<<endl;
}
