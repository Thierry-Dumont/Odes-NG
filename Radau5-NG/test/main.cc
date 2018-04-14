#define LOGRADAU5
#include <iostream>
#include <fstream>
#include "Radau5cc.hpp"
#include "AvcF.hpp"
#include "Oregonator.hpp"
#include "KPP.hpp"
#include "E5.hpp"
#include "BZ.hpp"
#include <time.h>
using namespace std;
int main()
{
  //---- define here the problem you want to treat:
  typedef Oregonator Fonc;
  //typedef KPP Fonc;
  //typedef AvcF Fonc;
  //typedef E5 Fonc;
  //typedef BZ Fonc;
  //---------------------------------------------- 
  cout.precision(17);

  static const int n=Fonc::n;
  double atol[n],rtol[n]; atol[0]=1.e-5; rtol[0]=1.e-5;

#ifdef ICC
  cout<<"Compiled with ICC."<<endl;
#endif
  cout<<n<<" equations."<<endl;
  
  Radau5cc<Fonc> Rad(true,atol,rtol);
  double h,t,xend;
  int nloops;

  double y[n];
  cout<<"initial time step?"; cin>>h;
  cout<<"integration time?";  cin>>xend;
  cout<<"how many loops on the problem?"; cin >>nloops;
  
  clock_t clkStart =  clock();// we measure execution time.
  for(int i=0;i<nloops;i++)//loop on the same problem.
    {
      Rad.rhs().init(y);
      Rad.setNmax(10000);
      t=0.0;
      try{
	Rad(h,t,xend,&(y[0]));
      }
      catch( OdesException )
	{
	  #ifdef LOGRADAU5
	  ofstream logfile; logfile.open("logfile");
	  logfile<<Rad.Log()<<endl;
	  logfile.close();
	  throw OdesException("abort! logfile created.");
	  #else
	  throw OdesException("abort!");
	  #endif
	  
	}
    }
  auto texec = static_cast<double>(clock() - clkStart)/CLOCKS_PER_SEC;
  
  cout<<endl<<endl;
  
  cout<<"nstep= "<<Rad.getNstep()<<"\nnjac= "<<Rad.getNJac()<<
    "\nnaccpt= "<<Rad.getNaccpt()<<"\nnreject= "<<Rad.getNrejct()<<
     "\nndec= "<<Rad.getNdec()<<endl;
  
  cout<<"first accepted step: "<<Rad.getfirstAcceptedStep()<<endl;
  
  cout<<"\nExecution time (milli seconds per loop): "<<1000.*texec/nloops<<endl;

#ifdef LOGRADAU5
  //std::cout<<Rad.Log()<<endl;// uncomment produces a lot of output.
  ofstream logfile; logfile.open("logfile");
  logfile<<Rad.Log()<<endl;
  logfile.close();
  cout<<"\nlogfile created."<<endl;
#endif
  cout<<endl<<"end."<<endl;
}
