#ifndef Extrap__h
#define Extrap__h
#include "OdesException.hpp"
#include "extrapcoeffs.hpp"
//////////////////////////////////////////////////////////////////////////////
/// Class (template of) for the prediction of values at the next step.
/// (see "Geometric ...." 2nd edition, pages 326-327).
///////////////////////////////////////////////////////////////////////////////
namespace odes{

  template<int n,int nsteps,class Double> struct Extrap
  {
    static const int n2=nsteps*nsteps;
    Double q[n2];
 
    //! constructor
    Extrap()
    {
      init();
    }
    //! destructor.
    ~Extrap(){}
    //! multiply coefficients by h.
    //! \param h
    void change_h(double h)
    {
      init();
      for(int i=0;i<n2;i++)
	q[i]*=h;
    }
    //! compute the predicted values.
    //! \param F  values of f()
    //! \param Y  unknows at substep
    //! \param u  values at preceding time step.
    void operator()(Double *F,Double *Y,Double *u)
    {
      for(int i=0;i<nsteps;i++)
	{
	  for(int j=0;j<n;j++)
	    Y[nsteps*i+j]=u[j];
	  for(int s=0;s<nsteps;s++)
	    for(int j=0;j<n;j++)
	      Y[nsteps*i+j]+=q[nsteps*s+j]*F[nsteps*s+j];
	}
	    
    }
private:
  void init()
  {
    // copy coefficients in q[].
    extrapcoeffs<nsteps,Double>(q);
  }
};


};
#endif
