#ifndef HESSENBERG__H
#define HESSENBERG__H
#include "fortranArray.hpp"
#include "fortranVector.hpp"
using namespace odes;
#define MIN(x,y) ((x)<(y)?(x):(y))
#define ABS(a) (((a) >= (0.0)) ? (a) : (-a))
/////////////////////////////////////////////////////////////////////////////
/// All what is necessary to compute with Hessenberg matrices (but not
/// the reduction of a full matrix to Hessenberg form (see Matrices.hpp).
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
/// Triangularization by Gaussian elimination of a Hessenberk Matrix with
/// lower Bandwidth ==1.
/// \brief Triangularization by Gaussian elimination of a Hessenberk Matrix.
/////////////////////////////////////////////////////////////////////////////
//! \param A the matrix (IN/OUT).
//! \param ip index of pivots.
//! \note return 0 ii everything went well, otherwise the rank wher the
//!        matrix A was found singular.
//! \note this is a C++ transcription of DECH (in H.&W. fortran code 
//!       decsol).
template<int n> int dech(fortranArray<n>& A,int ip[])
{
  int ret=0;
  ip[n-1]=1;
 
  if(n!=1)
    for(int k=1;k<n;k++)
      {
	int kp1=k+1,m=k;
	int na=MIN(n,1+k);// lb if fortran code; here lb==1.
	for(int i=kp1;i<=na;i++)
	  if(ABS(A(i,k))>ABS(A(m,k))) m=i;
	ip[k-1]=m;
	double t=A(m,k);
	if(m!=k)
	  {
	    ip[n-1]=-ip[n-1];
	    A(m,k)=A(k,k);
	    A(k,k)=t;
	  }
	if(t==0.0)
	  {
	    ret=k;
	    break;
	  }
	t=-1.0/t;
	for(int i=kp1;i<=na;i++)
	  A(i,k)*=t;
	for(int j=kp1;j<=n;j++)
	  {
	    t=A(m,j);
	    A(m,j)=A(k,j);
	    A(k,j)=t;
	    if(t!=0.0)
	      {
		for(int i=kp1;i<=na;i++)
		  A(i,j)+= A(i,k)*t;
	      }
	  }
      }
  if(A(n,n)==0.0) ret=n;
  return ret;
}

//////////////////////////////////////////////////////////////////////
/// Solution of a linear system where the matrix A was computed in
/// dech. 
//////////////////////////////////////////////////////////////////////
//! \param A the matrix computed in dech (IN)
//! \param ip the pivots obtained  in dech. (IN)
//! \param B the RHS on entry, the solution at the end (IN/OUT).
//! \note this is a transcription of the SOLH routine in H.& W. decsol. 
template<int n> void solh(const fortranArray<n>& A,const int ip[],
			 fortranVector& B)
{
  if(n!=1)
    {
      for(int k=1;k<n;k++)
	{
	  int kp1=k+1;
	  int M=ip[k-1];
	  double t=B(M);
	  B(M)=B(k);
	  B(k)=t;
	  int na=MIN(n,k+1);//LB==1.
	  for(int i=kp1;i<=na;i++)
	    B(i)+=A(i,k)*t;
	}
      for(int kb=1;kb<n;kb++)
	{
	  int km1=n-kb;
	  int k=km1+1;
	  B(k)/=A(k,k);
	  double t=-B(k);
	  for(int i=1;i<=km1;i++)
	    B(i)+=A(i,k)*t;
	}
    }
  B(1)/=A(1,1);
}

/////////////////////////////////////////////////////////////////////////////
/// Triangularization by Gaussian elimination of a Hessenberk Matrix with
/// lower Bandwidth ==1. The matrix, here is complex, but treated as a pair
/// of 2 real matrices.
/// \brief Triangularization by Gaussian elimination of a complex
/// \brief Hessenberg Matrix.
/// \note pure transciption of H. & W. fortran code.
/////////////////////////////////////////////////////////////////////////////
//! \param Ar the matrix,real part (IN/OUT).
//! \param Ai the matrix,imaginary part (IN/OUT).
//! \param ip index of pivots.
//! \note return 0 ii everything went well, otherwise the rank wher the
//!        matrix A was found singular.
//! \note this is a C++ transcription of DECHC (in H.&W. fortran code 
//!       decsol.
template<int n> int dechc(fortranArray<n>& Ar,fortranArray<n>& Ai,int ip[])
{
  int ret=0;
  ip[n-1]=1;
 
  if(n!=1)
    for(int k=1;k<n;k++)
      {
	int kp1=k+1,m=k;
	int na=MIN(n,1+k);// lb if fortran code; here lb==1.
	for(int i=kp1;i<=na;i++)
	  if(ABS(Ar(i,k))+ABS(Ai(i,k))>ABS(Ar(m,k))+ABS(Ai(m,k))) m=i;
	ip[k-1]=m;
	double tr=Ar(m,k),ti=Ai(m,k);
	if(m!=k)
	  {
	    ip[n-1]=-ip[n-1];
	    Ar(m,k)=Ar(k,k);
	    Ai(m,k)=Ai(k,k);
	    Ar(k,k)=tr; Ai(k,k)=ti;
	  }
	if(ABS(tr)+ABS(ti)==0.0)
	  {
	    ret=k;
	    break;
	  }
	double den=tr*tr+ti*ti;
	tr/=den;        ti/=-den;
	for(int i=kp1;i<=na;i++)
	  {
	    double prodr=Ar(i,k)*tr-Ai(i,k)*ti, prodi=Ai(i,k)*tr+Ar(i,k)*ti;
	    Ar(i,k)=-prodr;	                Ai(i,k)=-prodi;
	  }
	for(int j=kp1;j<=n;j++)
	  {
	    tr=Ar(m,j);        ti=Ai(m,j);
	    Ar(m,j)=Ar(k,j);   Ai(m,j)=Ai(k,j);
	    Ar(k,j)=tr;        Ai(k,j)=ti;
	    if(ABS(tr)+ABS(ti)!=0.0)
	      {
		if(ti==0.0)
		  for(int i=kp1;i<=na;i++)
		    {
		      double prodr=Ar(i,k)*tr,  prodi=Ai(i,k)*tr;
		      Ar(i,j)+=prodr;           Ai(i,j)+=prodi;
		    }
		else if(tr==0.0)
		  for(int i=kp1;i<=na;i++)
		    {
		      double prodr=-Ai(i,k)*ti,  prodi=Ar(i,k)*tr;
		      Ar(i,j)+=prodr;           Ai(i,j)+=prodi;
		    }
		else
		  for(int i=kp1;i<=na;i++)
		    {
		      double prodr= Ar(i,k)*tr-Ai(i,k)*ti;
		      double prodi= Ai(i,k)*tr+Ar(i,k)*ti;
		      Ar(i,j)+=prodr;           Ai(i,j)+=prodi;
		    }
	      }
	  }
      }
  if(ABS(Ar(n,n))+ABS(Ai(n,n))==0.0) ret=n;
  return ret;
}
//////////////////////////////////////////////////////////////////////
/// Solution of a linear system where the matrix A was computed in
/// dech. 
//////////////////////////////////////////////////////////////////////
//! \param Ar  matrix computed in dechc (IN)
//! \param Ai  matrix computed in dechc (IN)
//! \param ip the pivots obtained  in dech. (IN)
//! \param Br the RHS on entry (real part) the solution at the end (IN/OUT).
//! \param Bi the RHS on entry (imag. part) the solution at the end (IN/OUT).
//! \note this is a C++ transcription of SOLHC (in H.&W. fortran code 
//!       decsol)
template<int n> void solhc(const fortranArray<n>& Ar,const fortranArray<n>& Ai,
			   const int ip[],fortranVector& Br,fortranVector& Bi)
{
  if(n!=1)
    {
      for(int k=1;k<n;k++)
	{
	  int kp1=k+1;
	  int M=ip[k-1];
	  double tr=Br(M),    ti=Bi(M);
	  Br(M)=Br(k);        Bi(M)=Bi(k);
	  Br(k)=tr;           Bi(k)=ti;
	  int na=MIN(n,k+1);//LB==1.
	  for(int i=kp1;i<=na;i++)
	    {
	      double prodr=Ar(i,k)*tr-Ai(i,k)*ti;
	      double prodi=Ai(i,k)*tr+Ar(i,k)*ti;
	      Br(i)+=prodr;   Bi(i)+=prodi;
	    }
	}
      for(int kb=1;kb<n;kb++)
	{
	  int km1=n-kb;
	  int k=km1+1;
	  double den=Ar(k,k)*Ar(k,k) + Ai(k,k)*Ai(k,k);
	  double prodr= Br(k)*Ar(k,k)+Bi(k)*Ai(k,k);
	  double prodi= Bi(k)*Ar(k,k)-Br(k)*Ai(k,k);
	  Br(k)=prodr/den;    Bi(k)=prodi/den;
	  double tr= -Br(k);  double ti=-Bi(k);
	  for(int i=1;i<=km1;i++)
	    {
	      prodr=Ar(i,k)*tr-Ai(i,k)*ti;
	      prodi=Ai(i,k)*tr+Ar(i,k)*ti;
	      Br(i)+=prodr;   Bi(i)+=prodi;
	    }
	}
    }
  double den=Ar(1,1)*Ar(1,1) + Ai(1,1)*Ai(1,1);
  double prodr= Br(1)*Ar(1,1)+Bi(1)*Ai(1,1);
  double prodi= Bi(1)*Ar(1,1)-Br(1)*Ai(1,1);
  Br(1)=prodr/den;    Bi(1)=prodi/den;
}
////////////////////////////////////////////////////////////////////////////
/// Given a Matrix A, in Hessenberg form, as computed by DGEHRD (lapack)
/// with ilo=1 and ihi=n, and a vector X, compute:
///      Y=(Product of Householder reflexions) X.
///////////////////////////////////////////////////////////////////////////
//! \param A matrix computed in DGEHRD (with ilo=1 and ihi=n).
//! \param tau parameters computed in DGEHRD.
//! \param X (IN).
//! \param Y result.
//! \param direct if true apply Householder refelexions from 1 to n-2
//!               otherwise from n-2 down to 1.
//! \note that we use dgehrd (in Matrices.hpp) to reduce the Jacobian
//! matrix to Hessenberg form. We must adapt the computation to the
//! output of this routine.
template<int n> void householder(const fortranArray<n>& A,const double tau[],
				 const fortranVector& X,fortranVector& Y,
				 bool direct)
{
  Y=X;
  int jbeg,jend,jinc;
  if(direct)
    {
      jbeg=1;jend=n-1;jinc=1;
    }
  else
    {
      jbeg=n-2; jend=0;jinc=-1;
    }
  for(int j=jbeg;j!=jend;j+=jinc)
    {
      double s=Y(j+1);
      for(int i=j+2;i<=n;i++)
	s+=Y(i)*A(i,j);
      s*=tau[j-1];
      
      Y(j+1)-=s;
      for(int i=j+2;i<=n;i++)
	Y(i)-=s*A(i,j);
    }
}
#endif
