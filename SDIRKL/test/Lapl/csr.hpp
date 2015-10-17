#ifndef csr_h
#define csr_h
#include "Premat.hpp"
#include <iostream>
///////////////////////////////////////////////////////////////////////
/// A class of csr matrices.
//////////////////////////////////////////////////////////////////////
class CSRLineWise
{
  //to store a matrix in csr format
 protected:
  int order,nzero;
  int *ia,*ja;
  double *a;
  //
 
  typedef PrematLineWise PM;
  mutable int rhopass; mutable double the_rho;
 public:
  //! constructor
  CSRLineWise(){}
  //! init the matrix.
  //!  \param P a PrematLineWise matrix.
  void init(PrematLineWise& P);
  ~CSRLineWise();
  //! return integer ia array
  inline int* Ia() const {return ia;}
  //! return integer ja array
  inline int* Ja() const {return ja;}
  //! return array of double
  inline double* A() const {return a;}
  //! return "order" ie number of lines and columns.
  inline int Order() const {return order;}
  //! return number of non zero terms.
  inline int Nzeros() const {return nzero;}
  //! print the matrix.
  void print() const;

  void get_arrays(int* ix[],double* x[])
 {
   ix[0]=ia;ix[1]=ja;
   x[0]=a;
 }

};

#endif

