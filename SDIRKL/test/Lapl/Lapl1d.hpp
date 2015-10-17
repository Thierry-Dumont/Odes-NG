#ifndef Lapl1d__h
#define Lapl1d__h
#include "Premat.hpp"
#include "SuperLu.hpp"
#include "BuildMatt1d.hpp"
class Lapl1d
{
  SuperLu S;
  PrematLineWise P;
  int n;
public:
  Lapl1d(){}
  ~Lapl1d(){}
  void init(int _n, double _alpha)
  {
    n=_n;
    BuildMat(_n,P,_alpha);
    S.init(P);
    P.clear();
  }
  void apply(double RHS[],double X[])
  {
    S.solve(X,RHS);
  }

};
#endif
