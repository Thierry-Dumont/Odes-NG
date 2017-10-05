#ifndef fortranVector__h
#define fortranVector__h
#include "OdesException.hpp"
#include "fortranArray.hpp"
#include "AllocateDestroyVector.hpp"
#include <algorithm>
#include <string>
namespace odes
{
  /////////////////////////////////////////////////////////////////////////
  ///
  //// Vector class (of doubles). We want to use fortran indexing
  ///  (from  1 to...): this is just a wrapper around an array of doubles.
  ///  \brief Vector class (of doubles).
  ///
  /// Note: this is not a very *safe* class!.
  ///       (we are adult programmers :-) ).
  /////////////////////////////////////////////////////////////////////////
class fortranVector
{
#ifdef DEBUG
  std::string name;
#endif
  template<int k> friend class fortranArray;
  template<int k> friend class fortranVectorF;
protected:
  double *x;
  int size;
  bool deletable;
public:
  //! constructor
  //! \note we construct an empty class, not usable.
  fortranVector()
  {
    deletable=false; size=0; x=0;
  }
  //!  constructor
  //! \param k  size of the vector
  fortranVector(int k)
  {
    deletable=true;
    //x=new double[k];
    x=allocDoubleArray(k);
    size=k;
  }
  //! build a fortranVector.
  //! \param _x a vector of double
  //! \param _size size of the vector
  //! \note we "adopt" a vector, to be able to use fortran indexing;
  //! this is not safe!
  fortranVector(double *_x,int _size)
  {
    deletable=false;
    x=_x; size=_size;;
  }
  //! destructor
  ~fortranVector()
  {
#ifdef DEBUG
    if(x==0&&deletable) 
      throw OdesException("~fortranVector: not allocated");
#endif
    //if(deletable) delete[] x;
    if(deletable) destroyDoubleArray(x);
  }
  //! operator (): indexation 
  //! \param i index (between 1 and n (included)).
  inline double operator()(int i) const
  {
#ifdef DEBUG
    if(i<1||i>size)
      throw OdesException("fortranVector(), bad i",i,size);
#endif
    return *(x+i-1);
  }
  //! operator (): indexation 
  //! \param i index (between 1 and n (included)).
  //! \note returns a reference (can be used as lvalue).
  inline double& operator()(int i) 
  {
#ifdef DEBUG
    if(i<1||i>size)
      throw OdesException("fortranVector(&), bad i",i,size);
#endif
    return *(x+i-1);
  }  
  //! return a pointer to the vector of double.
  inline double* operator&() const {return x;}
  //! copy constructor
  //! \param V object to copy
  //! \note nothing is checked! not safe.
  inline void operator=(fortranVector V)
  {
    //std::copy(V.x,V.x+size,x);
    ASSUME_ALIGNED(x);ASSUME_ALIGNED(V.x);
#include "Ivdep.hpp"
    for(int i=0;i<size;i++)
      x[i]=V.x[i];
  }
  //! operator +=
  //! \param F to be added.
  inline void operator+=(fortranVector& F)
    {
      ASSUME_ALIGNED(x);ASSUME_ALIGNED(F.x);
#include "Ivdep.hpp"
      for(int i=0;i<size;i++)
	x[i]+=F.x[i];
    }
  //! set size
  //! _param  _size the size.
  inline void setsize(unsigned int _size){size=_size;}
  //! get the size
  inline int get_size() const {return size;}
  //! print
  void print(string s="") const
  {
    cout<<s<<" ";
    for(int i=0;i<size;i++)
      cout<<x[i]<<' ';
    cout<<endl;
  }
#ifdef DEBUG
  //! for debuging. Put a name.
  //! \param  _name the name.
  void setname(string _name){name=_name;}
#endif
};
  /////////////////////////////////////////////////////////////////////////
  ///
  //// Vector class (of doubles) OF FIXED SIZE n. 
  ///  We want to use fortran indexing
  ///  (from  1 to...): this is just a wrapper around an array of doubles.
  ///  \brief Vector class (of doubles).
  ///
  /// Note: this is not a very safe class!.
  /////////////////////////////////////////////////////////////////////////
  template<int n> class fortranVectorF: public fortranVector
  {
// #ifndef ON_STACK
//     double *y;
// #else
//     double y[n];
// #endif
    double *y;
    using fortranVector::x;
    using fortranVector::size;
    using fortranVector::deletable;
    
  public:
    //! constructor. Object created is usable.
    fortranVectorF(): fortranVector()
    {
// #ifndef ON_STACK
//       y=new double[n];
// #endif  
      y=allocDoubleArray(n);
      x=&(y[0]); size=n;deletable=false;
    }
    //! destructor
    ~fortranVectorF()
    {
      destroyDoubleArray(y);
// #ifndef ON_STACK
//       delete[] y;
// #endif
      //deletion is in fortranVector.
    }
    //! swith adress with an other  fortranVectorF<n>
    //! \param v we switch with v.
    inline void switchadr(fortranVectorF<n>& v)
      {
	double *xt,*yt;
	xt=v.x;yt=v.y;
	v.x=x;v.y=y;
	x=xt;y=yt;
      }
    //! put the same value everywhere.
    //! \param v the value.
    inline void operator=(double v)
    {
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=v;
    }
    //! copy an array of double.
    //! \param tab[] the array.
    //! \note we do not check the size of tab[].
    inline void operator=(double tab[])
    {
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=tab[i];
    }
    //! this<- X+Y.
    //!\param  X
    //!\param Y
    //! \note we do not make any check on X and Y
    inline void setsum2(fortranVector& X,fortranVector& Y)
    {
#ifdef DEBUG
      if(n!=X.get_size() ||n!=Y.get_size())
	throw OdesException("fortranVectorF: setsum2, sizes differ",
			       n,X.get_size(),Y.get_size());
#endif
      double *xx=&X,*yy=&Y;
      ASSUME_ALIGNED(xx);ASSUME_ALIGNED(yy);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
      	y[i]=xx[i]+yy[i];
    }
    //! operator=
    //! \param X *this=X.
    inline void operator=(fortranVectorF<n>& X)
    {
      ASSUME_ALIGNED(x);ASSUME_ALIGNED(X.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=X.y[i];
    }
    //! operator +=
    //! \param F to be added.
    inline void operator+=(fortranVectorF<n>& F)
    {
      ASSUME_ALIGNED(y);ASSUME_ALIGNED(F.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]+=F.y[i];
    }
    //! this = X +Y.
    //! \param X fortranVectorF<n>
    //! \param Y fortranVectorF<n>
    inline void sum(fortranVectorF<n>& X,fortranVectorF<n>& Y)
    {
      ASSUME_ALIGNED(y);ASSUME_ALIGNED(X.y);ASSUME_ALIGNED(Y.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]=X.y[i]+Y.y[i];
    }
    //! this = X +Y.
    //! \param X fortranVectorF<n>
    //! \param Y fortranVector<n>
    inline void sum(fortranVectorF<n>& X,fortranVector& Y)
    {
     ASSUME_ALIGNED(y);ASSUME_ALIGNED(X.y);ASSUME_ALIGNED(Y.x);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]=X.y[i]+Y.x[i];
    }
    //! this= a*Y
    //! \param a a double
    //! \param Y fortranVector<n>
    inline void a_Y(double a, const fortranVectorF<n>& Y)
    {
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]=a*Y.y[i];
    }
    //! this= X+ a1*Y1 + a2*Y2.
    //! \param X
    //! \param a1
    //! \param Y1
    //! \param a2
    //! \param Y2
    inline void x_lc2(fortranVectorF<n>& X, 
		 double a1, fortranVectorF<n>& Y1,
		 double a2, fortranVectorF<n>& Y2)
    {
      ASSUME_ALIGNED(y);ASSUME_ALIGNED(X.y);ASSUME_ALIGNED(Y1.y);
      ASSUME_ALIGNED(Y2.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]=X.y[i]+a1*Y1.y[i]+a2*Y2.y[i];

    }
    //! this = linear combination of 2 fortranVectorF<n>
    //! \param a coefficient
    //! \param x1 fortranVectorF<n>
    //! \param b coefficient
    //! \param x2 fortranVectorF<n>
    inline void lc2(double a,fortranVectorF<n>& x1,
		    double b,fortranVectorF<n>& x2)
    {
     ASSUME_ALIGNED(x);ASSUME_ALIGNED(x1.y);ASSUME_ALIGNED(x2.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=a*x1.y[i]+ b*x2.y[i];
    }
    //! this = linear combination of 3 fortranVectorF<n>
    //! \param a coefficient
    //! \param x1 fortranVectorF<n>
    //! \param b coefficient
    //! \param x2 fortranVectorF<n>
    //! \param c coefficient
    //! \param x3 fortranVectorF<n>
    inline void lc3(double a,fortranVectorF<n>& x1,
		    double b,fortranVectorF<n>& x2,
		    double c,fortranVectorF<n>& x3)
    {
      ASSUME_ALIGNED(x);ASSUME_ALIGNED(x1.y);ASSUME_ALIGNED(x2.y);
      ASSUME_ALIGNED(x3.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=a*x1.y[i]+ b*x2.y[i]+ c*x3.y[i];
    }
    //! this = linear combination of 4 fortranVectorF<n>
    //! \param a coefficient
    //! \param x1 fortranVectorF<n>
    //! \param b coefficient
    //! \param x2 fortranVectorF<n>
    //! \param c coefficient
    //! \param x3 fortranVectorF<n>
    //! \param d coefficient
    //! \param x4 fortranVectorF<n>
    inline void lc4(double a,fortranVectorF<n>& x1,
		    double b,fortranVectorF<n>& x2,
		    double c,fortranVectorF<n>& x3,
		    double d,fortranVectorF<n>& x4)
    {
     ASSUME_ALIGNED(x);ASSUME_ALIGNED(x1.y);ASSUME_ALIGNED(x2.y);
     ASSUME_ALIGNED(x3.y);ASSUME_ALIGNED(x4.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=a*x1.y[i]+ b*x2.y[i]+ c*x3.y[i]+d*x4.y[i];
    }
    //! this = linear combination of 5 fortranVectorF<n>
    //! \param a coefficient
    //! \param x1 fortranVectorF<n>
    //! \param b coefficient
    //! \param x2 fortranVectorF<n>
    //! \param c coefficient
    //! \param x3 fortranVectorF<n>
    //! \param d coefficient
    //! \param x4 fortranVectorF<n>
    //! \param e coefficient
    //! \param x5 fortranVectorF<n>
    inline void lc5(double a,fortranVectorF<n>& x1,
		    double b,fortranVectorF<n>& x2,
		    double c,fortranVectorF<n>& x3,
		    double d,fortranVectorF<n>& x4,
		    double e,fortranVectorF<n>& x5)
    {
      ASSUME_ALIGNED(x);ASSUME_ALIGNED(x1.y);ASSUME_ALIGNED(x2.y);
      ASSUME_ALIGNED(x3.y);ASSUME_ALIGNED(x4.y);ASSUME_ALIGNED(x5.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=a*x1.y[i]+ b*x2.y[i]+ c*x3.y[i]+d*x4.y[i]+e*x5.y[i];
    }
    //! this= X+ a1*Y1 + a2*Y2 + a3*Y3
    //! \param X
    //! \param a1
    //! \param Y1
    //! \param a2
    //! \param Y2
    //! \param a3
    //! \param Y3
    inline void x_lc3(fortranVectorF<n>& X, 
		      double a1, fortranVectorF<n>& Y1,
		      double a2, fortranVectorF<n>& Y2,
		      double a3, fortranVectorF<n>& Y3)
    {
      ASSUME_ALIGNED(x);ASSUME_ALIGNED(X.y);ASSUME_ALIGNED(Y1.y);
      ASSUME_ALIGNED(Y2.y);ASSUME_ALIGNED(Y3.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]=X.y[i]+a1*Y1.y[i]+a2*Y2.y[i]+a3*Y3.y[i];

    }
    //! this= X+ a1*Y1 + a2*Y2 + a3*Y3 + a4* Y4
    //! \param X
    //! \param a1
    //! \param Y1
    //! \param a2
    //! \param Y2
    //! \param a3
    //! \param Y3
    //! \param a4
    //! \param Y4    
    inline void x_lc4(fortranVectorF<n>& X, 
		      double a1, fortranVectorF<n>& Y1,
		      double a2, fortranVectorF<n>& Y2,
		      double a3, fortranVectorF<n>& Y3,
		      double a4, fortranVectorF<n>& Y4)
    {
      ASSUME_ALIGNED(x);ASSUME_ALIGNED(X.y);ASSUME_ALIGNED(Y1.y);
      ASSUME_ALIGNED(Y2.y);ASSUME_ALIGNED(Y3.y);ASSUME_ALIGNED(Y4.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]=X.y[i]+a1*Y1.y[i]+a2*Y2.y[i]+a3*Y3.y[i]+a4*Y4.y[i];

    }
    // this = a . fortranVectorF<n> + fortranVectorF<n>
    //! \param a coefficient
    //! \param X fortranVectorF<n>
    //! \param Y fortranVectorF<n>    
    inline void a_x_plus_y(double a,fortranVectorF<n>& X,
			   fortranVectorF<n>& Y)
    {
       ASSUME_ALIGNED(x);ASSUME_ALIGNED(X.y);ASSUME_ALIGNED(Y.y);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	x[i]=a*X.y[i]+Y.y[i];
    }

    //! copy in a vector of double.
    //! \param x[] the vector (result).
    //! \note we cannot check for the size of x[].
    inline void put(double _x[])
    {
      //std::copy(y,y+n,x);
      ASSUME_ALIGNED(x);ASSUME_ALIGNED(_x);
#include "Ivdep.hpp"
      for(int i=0;i<n;i++)
	y[i]=x[i];
    }
    //! print
    void print(string s="") const
    {
      cout<<s<<" ";
      for(int i=0;i<size;i++)
	cout<<y[i]<<' ';
      cout<<endl;
    }
  };
}
#endif
