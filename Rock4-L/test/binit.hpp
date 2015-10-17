#ifndef binit__h
#define binit__h
// initialize the RHS B (we solve du/dt=AU +B).
void binit(double b[],int n)
{
  for(int i=0;i<n;i++)
    b[i]=0;
}
#endif
