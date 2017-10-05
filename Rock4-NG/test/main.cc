#define LOGROCK4
#define ROCK4_OMP
#include <iostream>
#include <fstream>
#include "MacrosForCompilers.hpp"
#include "AllocateDestroyVector.hpp"
#include "Rock4.hpp"

#include "Lapl1OMP.hpp"
using namespace std;
using namespace odes;

int main()
{

  typedef Lapl1 F;

  static const int size=1500; 
  std::cout.precision(16);
  F Ff(size);

  Rock4<F> R(Ff);

  R.setTolerances(1.e-5,1.e-5);
  //
  double* y=allocDoubleArray(size);
  for(int rep=0;rep<1;rep++)
    {
      Ff.init(y);
  
      double t0=0.0, tend= 1., dt=tend;
      R(y,t0,tend,dt);
    }
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
#ifdef LOGROCK4
  R.Log().print();
#endif
  return 1;
}
