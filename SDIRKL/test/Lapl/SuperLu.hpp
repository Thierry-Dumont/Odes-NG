#ifndef SuperLu_h
#define SuperLu_h
#include "slu_ddefs.h"
#include "Premat.hpp"
class SuperLu
{
  typedef PrematLineWise Input;
  
  typedef pair<int,int> pos;
  typedef Input::map_type mat;
  int n,m,nnz;
  double *a,*rhs;
  int *rowind;
  int *colptr;
  int            *perm_c; /* column permutation vector */
  int            *perm_r; /* row permutations from partial pivoting */

  SuperMatrix A,Lf,Uf,B;
  superlu_options_t options;
  mem_usage_t    mem_usage;
  SuperLUStat_t stat;
  int info;
  int nrhs;
 public:
  SuperLu(){}
  void init(Input& M);
  SuperLu(Input& P){init(P);}
  ~SuperLu();
  void solve(double *B,double *sol)
    {
      for(int i=0;i<n;i++)
	B[i]=sol[i];
      solve(B);
    }
  void solve(double *BX);
};

#endif
