//#define LOGROCK4
#define ROCK2_OMP
#include <iostream>
#include <fstream>
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include "Rock2.hpp"
#include <time.h>
#include "FKPPOMP.hpp"
using namespace std;
using namespace odes;

int main()
{

  typedef FKPP F;

  int size;
  double h,xend;
  std::cout.precision(16);
  

  
  cout<<endl<<"Spatial discretization: how many points?"; cin>>size;
  cout<<"Initial time step?"; cin>>h;
  cout<<"Integration time?";  cin>>xend;
  cout<<endl<<"We integrate on t=[0,"<<xend<<"], with "<<
    size<<" points in x."<<endl<<endl;
  
  F Ff(size);
  Rock2<F> R(Ff);
  
  R.setTolerances(1.e-7,1.e-7);
  //
  double* y=allocDoubleArray(size);
  Ff.init(y);

   double t0=0.0;
  clock_t clkStart =  clock();// we measure execution time.
 
  try{
    R(y,t0,xend,h);
  }
  catch(odes::OdesException)
    {
      cout<<"If you got a message like 'Rock2: step too small',";
      cout<<" you probably need to restart with smaller initial time step";
      cout<<" and integration time."<<endl<<endl;;
#ifdef LOGROCK2
      cout<<R.Log()<<endl;
#endif
 
    }
  cout<<"ok, last time step: "<<R.LastAcceptedTimeStep()<<endl<<endl;
  
  auto texec = static_cast<double>(clock() - clkStart)/CLOCKS_PER_SEC;
 
  ofstream f; f.open("result");
  for(int i=0;i<size;i++)
    f<<y[i]<<endl;
  f.close();
  cout<<"Results at the end of the integration are in file './result'"<<
    endl<<endl;
  destroyDoubleArray(y);
  
#ifdef LOGROCK2
  cout<<R.Log()<<endl;
#endif  
  cout<<"Func called : "<<R.NbRhsComputed()<<endl;
  cout<<"Nb. stages  : "<<R.NbStages()<<endl;
  cout<<"Nb. steps   : "<<R.NbSteps()<<endl;
  cout<<"Nb. accepted: "<<R.NbAccepted()<<endl;
  cout<<"Nb. rejected: "<<R.NbRejected()<<endl;
  cout<<"Size        : "<<R.nbUnkn()<<endl;
  cout<<"Execution time (ms): "<<1000.*texec<<endl;

}
