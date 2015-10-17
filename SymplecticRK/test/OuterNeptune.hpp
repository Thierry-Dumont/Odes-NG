#ifndef OuterNeptune__h
#define OuterNeptune__h
#include <cmath>
template<class Double> class OuterNeptune
{
  inline Double dot(Double x[],Double y[],int k,int l) const
  {
    Double ret=0.0;
    for(int i=0;i<3;i++) ret+=x[3*k+i]*y[3*l+i];
    return ret;
  }
  inline Double diff2(Double x[],Double y[],int k,int l) const
  {
    Double ret=0.0;
    for(int i=0;i<3;i++)
      ret+=pow(x[3*k+i]-y[3*l+i],2);
    return ret;
  }
  inline Double norme_diff(Double x[],Double y[],int k,int l) const
  {
    return sqrt(diff2(x,y,k,l));
  }
 public:
  static const int nb=5;
  static const int n=6*nb;
  static const int nb1=nb-1;
  Double masses[nb],G;
  //constructor
  OuterNeptune()
  {
    //Jupiter:
    masses[0]=0.00095478610443;
    //Saturne:
    masses[1]=0.000285583733151; 
    //Uranus:
    masses[2]= 0.0000437273164546;
    //Neptune:
    masses[3]=0.0000517759138449;
     //Soleil:
    masses[4]=1.00000597682;
    //
    G=2.95912208286e-04;

  }
  inline void dqdt(Double p[],Double dq[]) const
  {
    for(int i=0;i<=nb1;i++)
      {
	Double mi=masses[i];
	for(int j=0;j<3;j++)
	  dq[3*i+j]=-p[3*i+j]/mi;
      }
  }
  inline void dpdt(Double q[],Double dp[]) const
  {
    for(int i=0;i<3*n;i++)
      dp[i]=0.0;
#ifdef INTELCOMP
#pragma parallel
#endif
    for(int l=1;l<=nb1;l++)
      for(int j=0;j<l;j++)
	{
	  Double c=G*masses[l]*masses[j]/pow(sqrt(diff2(q,q,l,j)),3);
	  for(int s=0;s<3;s++)
	    dp[3*l+s]+=c*(q[3*l+s]-q[3*j+s]);
	}
#ifdef INTELCOMP
#pragma parallel
#endif
    for(int l=0;l<nb1;l++)
      for(int j=l+1;j<=nb1;j++)
	{
	  Double c=-G*masses[j]*masses[l]/pow(sqrt(diff2(q,q,j,l)),3);
	  for(int s=0;s<3;s++)
	    dp[3*l+s]+=c*(q[3*j+s]-q[3*l+s]);
	}
  }
  inline void operator()(Double X[],Double Y[]) const
  {

    
    //symplify notations with references:
    Double *p=X+3*nb, *q=X;// q: position, p: momentum.
    Double *dp=Y+3*nb,*dq=Y;// q: position, p: momentum
    dqdt(p,dq);
    dpdt(q,dp);
    //

  }
  inline Double H(Double X[]) const //Hairer and Co. page 9..
  {
    //symplify notations:
    Double *p=X+3*nb,*q=X;// q: position, p: momentum.

    Double ret=0.0;
    for(int i=0;i<=nb1;i++)
      ret+=dot(p,p,i,i)/masses[i];
    ret/=2.0;
    Double s=0.0;
    for(int i=1;i<=nb1;i++)
      for(int j=0;j<i-1;j++)
	s+=masses[i]*masses[j]/norme_diff(q,q,i,j);
    ret-= G*s;
    return ret;
  }
  inline void init(Double U[])
  {
    Double *p=U+3*nb,*q=U;// q: position, p: momentum.
    //Jupiter:
    q[0]=-3.5023653; q[1]=-3.8169847; q[2]=-1.5507963;
    p[0]=0.00565429; p[1]=-0.00412490;p[2]=-0.00190589;
    //Saturne:
    q[3]=9.0755314; q[4]=-3.0458353; q[5]=-1.6483708;
    p[3]=0.00168318;  p[4]=0.00483525; p[5]=0.00192462;
    //Uranus:
    q[6]=8.3101420;q[7]=-16.2901086;q[8]=-7.2521278;
    p[6]=0.00354178;p[7]= 0.00137102; p[8]= 0.00055029;
    //Neptune:
    q[9]=11.4707666 ; q[10]=-25.7294829; q[11]= -10.8169456;
    p[9]=0.00288930 ; p[10]= 0.00114527; p[11]= 0.00039677;
    //Soleil:
    for(int i=0;i<3;i++)
      {q[12+i]=0.0; p[12+i]=0.0;}
    //momentum:
    for(int i=0;i<=nb1;i++)
      for(int j=0;j<3;j++)
	p[3*i+j]*=masses[i];
  }
};
#endif
