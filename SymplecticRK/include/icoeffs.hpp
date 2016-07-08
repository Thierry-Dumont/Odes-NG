#ifndef icoeffs__h
#define icoeffs__h
#include <initializer_list>
#include <type_traits>
using namespace std;
namespace odes{
  template<int nsteps,class Double> void icoeffs(Double* a,Double* b)
  {
    static_assert(nsteps>1 && nsteps<=16,"\n\n nsteps must be >1 and <=16.\n");
  }
  // template<> void icoeffs<2,double>(double* a,double* b)
  // {
  //   initializer_list<double> la={0.25,0.25-sqrt(3.)/6.0,
  // 				 0.25+sqrt(3.)/6.0,0.25};
  //   initializer_list<double> lb={0.5,0.5};
  //   copy(la.begin(),la.end(),a);
  //   copy(lb.begin(),lb.end(),b);
  // }
#include "../sage/generated_coeffs.hpp"
};
#endif
