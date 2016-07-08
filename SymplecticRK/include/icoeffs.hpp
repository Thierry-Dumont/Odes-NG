#ifndef icoeffs__h
#define icoeffs__h
#include <initializer_list>
using namespace std;
template<int nsteps,class Double> struct icoeffs
{
  void init(Double* a,Double* b)
  {
  }
};
template<> struct icoeffs<2,double>
{
  void init(double* a,double* b)
  {
    
    initializer_list<double> la={0.25,0.25-sqrt(3.)/6.0,
				 0.25+sqrt(3.)/6.0,0.25};
    initializer_list<double> lb={0.5,0.5};

    copy(la.begin(),la.end(),a);
    copy(lb.begin(),lb.end(),b);
  }
};
#endif
