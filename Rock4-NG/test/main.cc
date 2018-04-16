//#define LOGROCK4
#define ROCK4_OMP
#include <iostream>
#include <fstream>
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include "Rock4.hpp"
#include <time.h>
#include "FKPPOMP.hpp"
using namespace std;
using namespace odes;

int main()
{

  typedef FKPP F;

  static const int size=500; 
  std::cout.precision(16);
  F Ff(size);

  double h,xend;
  cout<<"initial time step?"; cin>>h;
  cout<<"integration time for each step?";  cin>>xend;
  Rock4<F> R(Ff);
  cout<<"We integrate on [0,"<<xend<<"]."<<endl<<endl;

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
      cout<<"If you got a message like 'Rock4: step too small',";
      cout<<" you probably need to restart with smaller initial time step";
      cout<<" and integration time."<<endl<<endl;;
#ifdef LOGROCK4
      cout<<R.Log()<<endl;
#endif
 
    }
 
  auto texec = static_cast<double>(clock() - clkStart)/CLOCKS_PER_SEC;
  cout<<"ok, last time step: "<<R.LastAcceptedTimeStep()<<endl<<endl;
 
  ofstream f; f.open("result");
  for(int i=0;i<size;i++)
    f<<y[i]<<endl;
  f.close();
  cout<<"Results at the end of the integration are in file './result'"<<
    endl<<endl;
  destroyDoubleArray(y);
#ifdef LOGROCK4
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
