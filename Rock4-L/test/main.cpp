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
  atol[0]=1.e-8; rtol[0]=1.e-8;
  double *y=allocDoubleArray(n);
  Radau5cc<FRadau5> Rad(true,atol,rtol);
  bool ok;

  cout<<endl<<"1) Find acceptable atol and rtol in Radau5, as small as possible:"<<endl;
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
	  cout<<"tol: "<<atol[0]<<" "<<rtol[0];//<<endl;
	}
    }
  while(!ok);

  
  cout<<endl<<"Ok, we have a reference solution at order 4, as precisely computed as possible with Radau5."<<endl<<endl;

  // Uncomment to print results.
  // cout<<"The solution, using Radau5: "<<endl;
  // print(y,n);
    
  cout<<"2) Now, we compute with Rock4, with decreasing time steps,";
  cout<<" and we compare the results with the 'best' solution obtained with Radau5:"<<endl<<endl;

  if(affine)
    cout<<"problem is AFFINE."<<endl;
  else
    cout<<"problem is LINEAR."<<endl;
  double *yR=allocDoubleArray(n);
   
  ofstream outR("OutR",ios::out|ios :: trunc);
  ofstream outF("OutF",ios::out|ios :: trunc);

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
      
 
     }
  // uncomment to print the results:
  // cout<<"n= "<<n<<endl;
  // for(int i=0;i<n;i++)
  //   outF<<yR[i]<<endl;
  // print(yR,n);
 
  destroyDoubleArray(yR);
  destroyDoubleArray(y);
  
  cout<<endl<<"The best solution computed with Rock4 is in ./OutF, error vs.time_step for Rock4 in ./OutR."<<endl;
  cout<<"You can just gnuplot these files (with logscale on both axis for ./OutR)."<<endl;
  
  cout<<endl;
}
