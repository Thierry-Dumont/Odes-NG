#include <iostream>
#include <fstream>
#include "SymplecticRK.hpp"
#include "OdesException.hpp"
#include "Outer.hpp"
#include "OuterUranus.hpp"
#include "OuterNeptune.hpp"
#include "Kepler2.hpp"
using namespace std;
using namespace odes;
int main()
{
  //typedef long double Double;
  typedef double Double;
  typedef Outer<Double> F;
  //typedef OuterUranus<Double> F;
  //typedef OuterNeptune<Double> F;
  //typedef Kepler2<Double> F;
  const int n=F::n;
  SymplecticRK<F,8,Double> Symp(100);
  Double u[n];
  Symp.rhs().init(u);
  //Symp.verify();
  Double h;
  cout.precision(53);
  //cout<<"h?"; cin>>h;
  ofstream fileout("result");
  cout<<"n= "<<n<<endl;
  h=1.0; //in days
  double years=200;
  bool firstep=true;
  for(int it=0;it<years*365.25*h;it++)
    {
      bool ok=Symp.step(h,u,firstep);
      if(!ok)
	{
	  cerr<<"last value: "<<Symp.lastdiff()<<endl;
	  throw OdesException("SymplecticRK::step: non convergence. h=",h);
	}
      for(int l=0;l<n;l++)
      	fileout<<u[l]<<' '; fileout<<'\n';
      
	fileout<<Symp.rhs().H(u)<<'\n';
      cout<<Symp.rhs().H(u)<<endl;
    }
  fileout.close();

  cout<<"end."<<endl;
}
