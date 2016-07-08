#include <iostream>
#include "protos_lapack.hpp"
#include "AllocateDestroyVector.hpp"
//#include "FRadau5.hpp"
//#include "Radau5cc.hpp"
#include "Rock4L.hpp"
#include "binit.hpp"
#include <fstream>
#include <cmath>
#include <string>
#include "classF.hpp"


using namespace std;
using namespace odes;

//! norm of the difference of 2 vectors.
double normdiff(int n,double x[],double y[])
{
  double ret=0.0;
  for(int i=0;i<n;i++)
    ret=pow((x[i]-y[i]),2);
  return sqrt(ret);
}
//! easy print a vector.
void print(double x[],int n)
{
  cout<<endl<<endl;
  for(int i=0;i<n;i++)
    cout<<x[i]<<" ";
  cout<<endl<<endl;
}
int main()
{
  //typedef FRadau5 F;
  static const int n=FRadau5::n;
  // Problem:
  static const bool affine=false;
  //
  double Tinteg=0.1;
  double *b=allocDoubleArray(n);//new double[n];
  //
  //Radau5 integration:
  //
  double atol[n],rtol[n]; 
  atol[0]=1.e-25; rtol[0]=1.e-25;
  double *y=allocDoubleArray(n);
  Radau5cc<FRadau5> Rad(true,atol,rtol);
  bool ok;
  do
    {
      // find acceptable atol and rtol, as small as possible.
      try
	{
	 
	  ok=true;
	  Rad.setTestPolicy(true,atol,rtol);
	  double t=0.0,h=Tinteg;
	  Rad.rhs().init(y);
	  Rad(h,t,Tinteg,&(y[0]));
	}
      catch(GenericException)
	{
	  ok=false;
	  atol[0]*=2; rtol[0]*=2;
	  cout<<"tol: "<<atol[0]<<" "<<rtol[0]<<endl;
	}
    }
  while(!ok);
  //ok, we have a solution at order 4, as precise as possible with Radau5.
  print(y,n);

  
  //
  ofstream f("gplot",ios::out|ios :: trunc);
  ofstream outR("OutR",ios::out|ios :: trunc);
  ofstream outF("OutF",ios::out|ios :: trunc);

  if(affine)
    cout<<"AFFINE"<<endl;
  else
    cout<<"LINEAR"<<endl;
  double *yR=allocDoubleArray(n);

  //compute using Rock, with decreasing time steps:
  for(int l=2;l<10;l++)
    {
      int np=pow(2,l);
      double dt=Tinteg/(double) np;
      double t0=0;

      F monF(n);
      Rock4L<F,affine> MonRock(monF);
     
      Rad.rhs().init(yR);
      cout<<"yR: ";
      for(int i=0;i<n;i++)
	cout<<yR[i]<<" ";cout<<endl;
      binit(b,n);
      if(affine)
	MonRock(yR,b,t0,Tinteg,dt);
      else
	MonRock(yR,t0,Tinteg,dt);


      cout<<endl<<"stages: "<<MonRock.NbStages()<<"  Nb. Pas: "<<
	MonRock.NbSteps()<<endl;
   
      //cout<<"t0,tend,dt "<<t0<<" "<<Tinteg<<" "<<dt<<endl;

      // difference between Radau5 computation and Rock4 computation:
      cout<<"ROCK err: "<<normdiff(n,yR,y)<<" dt: "<<dt<<endl;
 
      outR<<log10(dt)<<" "<<log10(normdiff(n,yR,y))<<endl;
      
      cout<<"n= "<<n<<endl;
      for(int i=0;i<n;i++)
	outF<<yR[i]<<endl;;
     }
  print(yR,n);
  f.close();
 
  destroyDoubleArray(yR);
  destroyDoubleArray(y);
  cout<<"end"<<endl;
}
