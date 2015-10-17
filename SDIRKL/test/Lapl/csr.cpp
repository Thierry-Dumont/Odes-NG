#include "csr.hpp"

#include <iostream>
using namespace std;
void CSRLineWise::init(PrematLineWise& P)
{
  order=P.max_indice_lig();
  nzero=P.ncoeffs();
  ia=new int[order+1];
  ja=new int[nzero];
  a=new double[nzero];
 
  PM::mat::const_iterator i=P.Mmap().begin();
  if(P.first(i)!=1)
    {
      cout<<P.first(i)<<endl;
      throw "not 1";
    }

  ia[0]=1; 
  
  int pja=0,curia=1; 
  for(;i!=P.Mmap().end();i++)
    {
      int cc=P.first(i);
      if(cc!=curia)
	{
	  for(int l=curia+1;l<=cc;l++)
	    ia[l-1]=pja+1;
	  curia=cc;
	}
      ja[pja]=P.second(i);
      a[pja++]=i->second;
    }
  ia[curia]=pja+1;
  rhopass=0;
}
CSRLineWise::~CSRLineWise()
{ 
  delete[] ia; delete[] ja; delete[] a;

}
void CSRLineWise::print() const
{
  cout<<"order= "<<order<<" nzeros= "<<nzero<<endl;
  for(int i=1;i<=order;i++)
    {
      cout<<"line: "<<i<<" start at: "<<ia[i-1]<<" end at:"<<ia[i]-1<<endl;
      for(int j=ia[i-1];j<=ia[i]-1;j++)
	cout<<'('<<ja[j-1]<<' '<<a[j-1]<<')'<<' ';
      cout<<endl;
    }
}



