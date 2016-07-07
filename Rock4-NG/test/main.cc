//#define USE_TBB
//#define Rock4_history
#include <iostream>
#include <fstream>
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include "Rock4.hpp"
#ifdef USE_TBB
#include "ToolTbb.hpp"
#include "LaplTBB.hpp"
#endif
//#include "ToolOMP.hpp"
#include "Lapl1OMP.hpp"
using namespace std;
using namespace odes;

int main()
{
#ifdef USE_TBB
  tbb::task_scheduler_init init;
#endif
#ifdef USE_TBB
  typedef LaplTBB F;
#else
  typedef Lapl1 F;
#endif
  static const int size=1500; 
  std::cout.precision(16);
  F Ff(size);
#ifdef USE_TBB
  Rock4<F,ToolTbb<F> > R(Ff);
#else
  Rock4<F> R(Ff);
#endif
  R.setTolerances(1.e30,1.e30);
  //
  double* y=allocDoubleArray(size);
  for(int i=0;i<size/2;i++)
    //y[i]=1.+(double) (i+1)/(double) size;
    y[i]=1.;
  for(int i=size/2;i<size;i++) y[i]=0.;
  //   cout<<y[i]<<" "; cout<<endl;
  double t0=0.0, tend= 1., dt=tend;
  R(y,t0,tend,dt);
  cout<<"ok, last time step: "<<R.LastAcceptedTimeStep()<<endl;
 
  ofstream f; f.open("result");
  for(int i=0;i<size;i++)
    f<<y[i]<<endl;
  f.close();
  destroyDoubleArray(y);
  cout<<"Func called : "<<R.NbRhsComputed()<<endl;
  cout<<"Nb. stages  : "<<R.NbStages()<<endl;
  cout<<"Nb. steps   : "<<R.NbSteps()<<endl;
  cout<<"Nb. accepted: "<<R.NbAccepted()<<endl;
  cout<<"Nb. rejected: "<<R.NbRejected()<<endl;
  cout<<"Size        : "<<R.nbUnkn()<<endl;
#ifdef Rock4_history
  R.get_history().print();
  cout<<endl<<R.get_history().lastSuccess(0)<<endl;
#endif
  return 1;
}
