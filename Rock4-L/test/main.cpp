#include <iostream>
#include "protos_lapack.hpp"
#include "AllocateDestroyVector.hpp"
#include "FRadau5.hpp"
#include "Radau5cc.hpp"
#include "Rock4L.hpp"
#include "binit.hpp"
#include <fstream>
#include <cmath>
#include <string>
#include "F.hpp"


using namespace std;
using namespace odes;
//
// We compare solutions with Radau5 and solutions with Rock4L.
//
//! norm of the difference of 2 vectors.
double normdiff(int n,double x[],double y[])
{
  double ret=0.0;
  for(int i=0;i<n;i++)
    ret=pow((x[i]-y[i]),2);
  return sqrt(ret);
}
//! easy print of a vector.
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
  double Tinteg=0.5;
  double *b=allocDoubleArray(n);//new double[n];
  //
  //Radau5 integration:
  //
  double atol[n],rtol[n]; 
  atol[0]=1.e-10; rtol[0]=1.e-10;
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
      catch(OdesException)
	{
	  ok=false;
	  atol[0]*=2; rtol[0]*=2;
	  cout<<"tol: "<<atol[0]<<" "<<rtol[0]<<endl;
	}
    }
  while(!ok);
  //ok, we have a solution at order 4, as precise as possible with Radau5.
  cout<<"Solution, using Radau5: "<<endl;
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

  //compute using Rock, with decreasing time steps, and compare with the
  // "best" solution obtained by Radau5:

  double dt = Tinteg/2.;
  for(int l=0;l<30;l++)
    {
      dt/= 1.25;
      double t0=0;

      F monF(n);
      Rock4L<F,affine> MonRock(monF);
     
      Rad.rhs().init(yR);
      binit(b,n);
      if(affine)
	MonRock(yR,b,t0,Tinteg,dt);
      else
	MonRock(yR,t0,Tinteg,dt);


      cout<<endl<<"stages: "<<MonRock.NbStages()<<"  Nb. steps: "<<
	MonRock.NbSteps()<<endl;
   
       // difference between Radau5 computation and Rock4 computation:
      cout<<"ROCK err: "<<normdiff(n,yR,y)<<" dt: "<<dt<<endl;
 
      outR<<dt<<" "<<normdiff(n,yR,y)<<endl;
      
      cout<<"n= "<<n<<endl;
      for(int i=0;i<n;i++)
	outF<<yR[i]<<endl;;
     }
  print(yR,n);
  f.close();
 
  destroyDoubleArray(yR);
  destroyDoubleArray(y);
  cout<<"The solution is in OutF, error/time_step in OutR."<<endl;
  cout<<"end"<<endl;
}
