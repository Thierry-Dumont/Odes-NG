void BuildMat(int n,PrematLineWise& P,double gamma)
{
 
  P.clear();
  double h1=1./(double) (n-1),h2=h1*h1,uh2=gamma/h2,duh2=2*uh2;
 
  for(int i=2;i<=n-1;i++)
    {
      P(i,i+1)=-uh2;
      P(i,i-1)=-uh2;
      P(i,i)=   duh2+1.;
	
    }
  P(1,1)= duh2+1.;
  P(1,2)=-duh2;

  P(n,n)=   duh2+1.;
  P(n,n-1)=-duh2;
  

  
}
