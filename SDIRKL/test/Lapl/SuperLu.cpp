
#include "SuperLu.hpp"
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include "GenericException.hpp"

void SuperLu::init(Input& Mm)
{
  
  n= Mm.max_indice_lig();
  m= n;
  nnz=Mm.ncoeffs();
  
  a=new double[Mm.ncoeffs()];

  rowind= new int[ Mm.ncoeffs()];

  colptr=new int[n+1];
  colptr[n]= Mm.ncoeffs();
  
  double *ap=a; int *rowi=rowind; 
  Input::mat mt=Mm.Mmap();
  Input::mat::const_iterator i=mt.begin();
  int l=0, colc=-1;

  for(;i!=mt.end();i++)
    {
      if(i->first.second-1!=colc)
	colptr[++colc]=l;
      *rowi++=i->first.first-1;
      *ap++=i->second;
      ++l;
    }

  set_default_options(&options);

  /* Now we modify the default options to use the symmetric mode. */
  options.SymmetricMode = YES;
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = 0.001;
  //------------------
  nrhs=1;
  set_default_options(&options);
  dCreate_CompCol_Matrix(&A, m, n,  Mm.ncoeffs(), a, rowind, colptr,
			 SLU_NC, SLU_D, SLU_GE);
//   dCreate_CompCol_Matrix(&A, m, n,  Mm.ncoeffs(), a, rowind, colptr,
// 			 SLU_NR, SLU_D, SLU_GE);

  //------------------
  dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
  perm_c=new int[n];
  perm_r=new int[n];

  /* Initialize the statistics variables. */
  StatInit(&stat);
  // factorization 
  B.ncol = 0; 


 /* cout<<"TEST = lda : "<<Bstore->lda<<" et row : "<< A.nrow <<endl;

  cout<<"       A.nrow : " <<A.nrow<<endl;
  cout<<"      B->Stype : " <<B.Stype<<" et SLU_DN = "<<SLU_DN<<endl;
  cout<<"      B->Dtype  : " <<B.Dtype<<" et SLU_D = "<< SLU_D<<endl;
  cout<<"      B->Mtype  : " <<B.Mtype<<" et SLU_GE = "<< SLU_GE<<endl;
  info=0;*/
  dgssv(&options, &A, perm_c, perm_r, &Lf, &Uf, &B, &stat, &info);

  // cout<<"Info = : "<<info<<endl;
 

    if ( info == 0 ) 
      {
	//	NCformat *Ustore;SCformat *Lstore;
	/* This is how you could access the solution matrix. */
        //double *sol = (double*) ((DNformat*) B.Store)->nzval; 
	
	 /* Compute the infinity norm of the error. */
	//dinf_norm_error(nrhs, &B, xact);

	//	std::cout<<"-SuperLu factorization:"<<std::endl;
	//	Lstore = (SCformat *) Lf.Store;
	//	Ustore = (NCformat *) Uf.Store;
	// 	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
	//	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
	//	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
	
	//dQuerySpace(&Lf, &Uf, &mem_usage);
	//	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	//      mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	//      mem_usage.expansions);
	delete[] a; delete[] rowind; delete[] colptr;
      }
    else
      throw GenericException("SuperLu::SuperLu, info = ",info);

}

//! destructor.
SuperLu::~SuperLu()
{
  delete[] perm_c;delete[] perm_r;
  Destroy_SuperMatrix_Store(&B);
 }


void SuperLu::solve(double *x)
{
  //double *x=xx.data();
  dCreate_Dense_Matrix(&B, m, nrhs, x, m, SLU_DN, SLU_D, SLU_GE);
  B.ncol =1;
  dgstrs(NOTRANS,&Lf,&Uf, perm_c, perm_r,&B,&stat,&info);
  if ( info != 0 )
    throw GenericException("SuperLu::solve, dgtrs: info = ",info);
}

