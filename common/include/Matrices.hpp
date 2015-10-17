#ifndef Matrices__h
#define Matrices__h
#include "fortranArray.hpp"
#include "fortranComplexArray.hpp"
#include "fortranRectangularArray.hpp"
#include "fortranRectangularComplexArray.hpp"
#include "fortranVector.hpp"
#include "AllocateDestroyVector.hpp"
#include "Hessenberg.hpp"
#include "GenericException.hpp"
//#include <cblas.h>
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
/////////////////////////////////////////////////////////////////////////////
/// Compute all matrix operations:
///       -compute Jacobian matrix
///       -....... Factorisations.
///       - solve systems.
/// Template for full matrix, with  specializations for banded matrices
/// and full matrix, transformed in Hessenberg form.
/// Trying to instanciate a banded matrix matrix in Hessenberg form
/// is impossible and should be capturated with "compat.hpp"
/// \brief Compute all matrix operations.
///////////////////////////////////////////////////////////////////////////
namespace odes
{
  //! full matrices.
  //case full== true (ie: n-nsub==1 and n-nsup==1), and Hessenberg=false. 
  template<bool full,bool Hessengerg,int n,int nsub,int nsup> class Matrices
  {
  public:
    typedef fortranArray<n> MatrixReal;
    typedef fortranComplexArray<n> MatrixComplex;
   
  private:
    int ipivr[n],ipivc[n];
    MatrixReal E1;
    MatrixComplex E2R;
    fortranVectorF<2*n> Z2N;// do not re-use Z2 in slvrad like in Hairer.
  protected:
    bool calhes;
    //§ constructor
    Matrices(){}
    //! build and factorize "real" matrix
    //! \param fac1 : we add fac1*I to the Jacobian.
    //! \param Jac the jacobian.
    inline void decomr(double fac1,MatrixReal& Jac)
    {
      E1.equal_minus(Jac);
      E1.addDiag(fac1);
      int nn=n,info;
      dgetrf_(&nn,&nn,&E1,&nn,&(ipivr[0]),&info);
      if(info!=0)
	throw GenericException("odes::Matrices::decomr dgetrf,info=",info);
    }
    //! build and factorize "complex" matrix
    //! \param alpha : we add alpha*I to the real part of the Jacobian.
    //! \param beta  : we add alpha*I to the imaginary part of the Jacobian.
    //! \param Jac the jacobian.
    inline void decomc(double alpha,double beta,const MatrixReal& Jac)
    {

      for(int j=1;j<=n;j++)
#include "Ivdep.hpp"
	for(int i=1;i<=n;i++)
	  E2R.set(i,j,-Jac(i,j));
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)
	{
	  E2R.Re(i,i)+=alpha;
	  E2R.Im(i,i)=beta;
	}
      int nn=n,info;
      zgetrf_(&nn,&nn,&E2R,&nn,&(ipivc[0]),&info);
      if(info!=0)
	throw GenericException("odes::Matrices::decomc zgetrf,info=",info);
      
    }
    //! solve "real" part.
    //! \param Z IN/OUT: 2nd member (IN); result (OUT).
    //! \param  Jac the Jacobian matrix.
    inline void solvereal(fortranVectorF<n>& Z,const MatrixReal& Jac)
    {
      int nn=n,un=1,ier;
      char notrans='n';
      dgetrs_(&notrans,&nn,&un,&E1,&nn,&(ipivr[0]),&Z,&nn,&ier);
      if(ier!=0)
	throw GenericException("odes::Matrices::solvereal, dgetrs,ier=",ier);
    }
    //! solve "imaginary" part.
    //! \param Zr IN/OUT: RHS, real part (IN); result,real part (OUT).  
    //! \param Zi IN/OUT: RHS, imag. part (IN); result,imag part (OUT).
    //! \param  Jac the Jacobian matrix.
    inline void solvecomplex(fortranVectorF<n>& Zr,fortranVectorF<n>& Zi,
    			     const fortranArray<n>& Jac)
    {
#include "Ivdep.hpp"
      for(int i=n;i>=1;i--)
    	{
    	  Z2N(2*i-1)=Zr(i);
    	  Z2N(2*i)=Zi(i);
    	}	  
      int nn=n,un=1,ier; char notrans='n';
      zgetrs_(&notrans,&nn,&un,&E2R,&nn,&(ipivc[0]),&Z2N,&nn,&ier);
      if(ier!=0)
    	throw GenericException("odes::Matrices::solvecomplex, zgetrs,ier=",ier);
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)
    	{
    	  Zi(i)=Z2N(2*i);
    	  Zr(i)=Z2N(2*i-1);
    	}
    }
    //! begining of line i.
    inline int Ibegin(int i) const {return 1;}
    //! end of line i.
    inline int Iend(int i) const {return n;}
  };
  ///////////////////////////////////////////////////////////////////////////
  //! banded matrices, with nsub subdiagonals (nsub>=0) and nsup (>=0)
  //! super diagonals.
  //! storage is adapted to lapack  banded matrices routines.
  //case full== false (ie: n-nsub>1 or n-nsup>1) and Hessenberg=false.
  /////////////////////////////////////////////////////////////////////////
  template<int n,int nsub,int nsup> class Matrices<false,false,n,nsub,nsup>
  {
  public:
    typedef fortranRectangularArray<n,nsub,nsup> MatrixReal;
    typedef fortranRectangularComplexArray<n,nsub,nsup> MatrixComplex;

  private:
    int ipivr[n],ipivc[n];
    MatrixReal E1;
    MatrixComplex E2R;
    fortranVectorF<2*n> Z2N;// do not re-use Z2 in slvrad like in Hairer.
    static const int ldab=2*nsub+nsup+1;
  protected:
    bool calhes;
    //! constructor.
    Matrices()
    {
      if(nsub>=n && nsup>=n)
	throw GenericException("Matrices (banded): incorrect nsub or nsup", 
			       "nsub=",nsub,"nsup=",nsup,"n=",n);
    }
    //! see full matrix case.
    inline void decomr(double fac1,const MatrixReal& Jac)
    {
      E1.equal_minus(Jac);
      E1.addDiag(fac1);
      int nn=n,knsub=nsub,knsup=nsup,lldab=ldab,info;
      dgbtrf_(&nn,&nn,&knsub,&knsup,&E1,&lldab,&(ipivr[0]),&info);
      if(info!=0)
	throw GenericException("odes::Matrices::decomr dgbtrf,info=",info);
    }
    //! see full matrix case.
    inline void decomc(double alpha,double beta,const MatrixReal& Jac)
    {
      for(int j=1;j<=n;j++)
	//#include "Ivdep.hpp" //gcc does not like ivdep here.
	for(int i=MAX(1,j-nsup);i<=MIN(n,j+nsub);i++)
	  E2R.set(i,j,-Jac(i,j));
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)
	{
	  E2R.Re(i,i)+=alpha;
	  E2R.Im(i,i)=beta;
	}
      int nn=n,knsub=nsub,knsup=nsup,lldab=ldab,info;
      zgbtrf_(&nn,&nn,&knsub,&knsup,&E2R,&lldab,&(ipivc[0]),&info);
      if(info!=0)
	throw GenericException("odes::Matrices::decomc zgetrf,info=",info);
    }
    //! see full matrix case.
    inline void solvereal(fortranVectorF<n>& Z,const MatrixReal& Jac)
    {
      int nn=n,knsub=nsub,knsup=nsup,lldab=ldab,un=1,ier;
      char notrans='n';
      dgbtrs_(&notrans,&nn,&knsub,&knsup,&un,&E1,&lldab,&(ipivr[0]),
	      &Z,&nn,&ier);
      if(ier!=0)
	throw GenericException("odes::Matrices::slvrad, dgetrs,ier=",ier);
    }
    //! see full matrix case.
    inline void solvecomplex(fortranVectorF<n>& Zr,fortranVectorF<n>& Zi,
			     const MatrixReal& Jac)
    {
#include "Ivdep.hpp"
      for(int i=n;i>=1;i--)
	{
	  Z2N(2*i-1)=Zr(i);
	  Z2N(2*i)=Zi(i);
	}
      int nn=n,knsub=nsub,knsup=nsup,lldab=ldab,un=1,ier; char notrans='n';
      zgbtrs_(&notrans,&nn,&knsub,&knsup,&un,&E2R,&lldab,&(ipivc[0]),
	      &Z2N,&nn,&ier);
      if(ier!=0)
	throw GenericException("odes::Matrices::slvrad, zgetrs,ier=",ier);
#include "Ivdep.hpp"
      for(int i=1;i<=n;i++)
	{
	  Zi(i)=Z2N(2*i);
	  Zr(i)=Z2N(2*i-1);
	}
    }
    //! begining of line i.
    inline int Ibegin(int i) const {return MAX(1,i-nsup);}
    //! end of line i.
    inline int Iend(int i) const {return MIN(n,i+nsub);}
  };
  ////////////////////////////////////////////////////////////////////////////
  //! Full matrices, Hessenberg=true.
  ///////////////////////////////////////////////////////////////////////////
  template<int n,int nsub,int nsup> class  Matrices<true,true,n,nsub,nsup>
  {
  public:
    typedef fortranArray<n> MatrixReal;
    typedef fortranComplexArray<n> MatrixComplex;

  private:
    double *work; int lwork;
    int ipr[n],ipi[n]; double tau[n];
    MatrixReal E1;
    MatrixReal E2R,E2I;
    fortranVectorF<n> Cr,Ci,Fr,Fi;
  protected:
    bool calhes;
    //! constructor.
    Matrices()
    {
      //work=new double[n];
      work=allocDoubleArray(n);
      lwork=n;
    }
    ~Matrices()
    {
      //delete[] work;
      destroyDoubleArray(work);
    }
    inline void decomr(double fac1,MatrixReal& Jac)
    {
      if(calhes)
	{
	  int ilo=1,ihi=n,lda=n,info;
	  int k=-1;
	  //find the optimum size for array work/
	  int nn=n;
	  dgehrd_(&nn,&ilo,&ihi,&Jac,&lda,tau,work,&k,&info);
	 
	  k=work[0];
	  if(lwork<k)
	    {
	      if(work!=0) delete[] work;
	      lwork=k;
	      work=new double[lwork];
	    }
	  dgehrd_(&nn,&ilo,&ihi,&Jac,&lda,tau,work,&k,&info);
	  if(info!=0)
	    throw GenericException("odes::Matrices::Matrices, decomr: info=",
				   info);
	  calhes=false;
	}
      E1.equal_minus(Jac);
      E1.addDiag(fac1);
      int info=dech<n>(E1,ipr);
      if(info!=0)
	throw GenericException("odes::Matrices::Matrices, decomr (dech): ",
			       info);
  
    }
    inline void solvereal(fortranVectorF<n>& Z,const MatrixReal& Jac)
    {
      householder<n>(Jac,tau,Z,Cr,true);
      solh<n>(E1,ipr,Cr);
      householder<n>(Jac,tau,Cr,Z,false);
    }
    inline void decomc(double alpha,double beta,const MatrixReal& Jac)
    {
      // here, Jac is in Hessenberg form.
#include "Ivdep.hpp"
      for(int j=1;j<n;j++)
      	{
      	  int j1=j+1;
      	  E2R(j1,j)=-Jac(j1,j);
      	  E2I(j1,j)=0.0;
      	}
      for(int j=1;j<=n;j++)
	{
#include "Ivdep.hpp"
	  for(int i=1;i<=j;i++)
	    {
	      E2I(i,j)=0.0; E2R(i,j)=-Jac(i,j);
	    }
	  E2R(j,j)+=alpha; E2I(j,j)=beta;
	}
      dechc(E2R,E2I,ipi);
    }
    inline void solvecomplex(fortranVectorF<n>& Zr,fortranVectorF<n>& Zi,
			     const fortranArray<n>& Jac)
    {
      householder<n>(Jac,tau,Zr,Cr,true);
      householder<n>(Jac,tau,Zi,Ci,true);
      solhc<n>(E2R,E2I,ipi,Cr,Ci);
      householder<n>(Jac,tau,Cr,Zr,false);
      householder<n>(Jac,tau,Ci,Zi,false); 
    }
    inline int Ibegin(int i) const {return 1;}
    inline int Iend(int i) const {return n;}
    };
    //---------------------------------------------------------------------
    //! Banded matrices, Hessenberg=true: error.
  template<int n,int nsub,int nsup> class Matrices<false,true,n,nsub,nsup>
  {
  public:
    typedef fortranArray<n> MatrixReal;
    typedef fortranComplexArray<n> MatrixComplex;
 
  protected:
    bool calhes;
    //! constructor.
    Matrices()
    {
      throw GenericException("Matrices cannot use Hessenberg option with",
  			     "banded Matrices");
    }
    inline void decomr(double fac1,const MatrixReal& Jac){}
    inline void decomc(double alpha,double beta,const MatrixReal& Jac){}
    inline void solvereal(fortranVectorF<n>& Z,const MatrixReal& Jac){}
    inline void solvecomplex(fortranVectorF<n>& Zr,fortranVectorF<n>& Zi,
  			     const fortranArray<n>& Jac){}
    inline int Ibegin(int i) const {return 1;}
    inline int Iend(int i) const {return n;}
  };
}
#endif
