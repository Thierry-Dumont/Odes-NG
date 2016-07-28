#ifndef Extrap__h
#define Extrap__h
#include "OdesException.hpp"
#include "extrapcoeffs.hpp"
//////////////////////////////////////////////////////////////////////////////
/// Class (template of) for the prediction of values at the next step.
/// (see "Geometric ...." 2nd edition, pages 326-327).
/// This is the method described at page 326 of Hairer and Co., based on an
/// extrapolation based on the collocation polynomial.
///////////////////////////////////////////////////////////////////////////////
namespace odes{

  template<int n,int nsteps,class Double> struct Extrap
  {
    static const int n2=nsteps*nsteps;
    Double q[n2];
 
    //! constructor
    Extrap()
    {
      // copy coefficients in q[].
      extrapcoeffs<nsteps,Double>(q);
    }
    //! destructor.
    ~Extrap(){}
    //! multiply coefficients by h.
    //! \param h
    void change_h(double h)
    {
       // copy coefficients in q[].
      extrapcoeffs<nsteps,Double>(q);
      //
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
	  for(int l=0;l<n;l++)
	    Y[nsteps*i+l]=u[l];
	  for(int j=0;j<nsteps;j++)
	    {
	      double beta=q[nsteps*i+j];
	      for(int l=0;l<n;l++)
		Y[nsteps*i+l]+=beta*F[nsteps*j+l];
	    }
	}
	    
    }

};


};
#endif
