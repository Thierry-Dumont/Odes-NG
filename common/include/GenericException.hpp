#ifndef GenericException_h
#define GenericException_h
#include <iostream>
using namespace std;
/////////////////////////////////////////////////////////////////////////////
///
/// A simple exception class.
/// Instantiate it with 1 to 6 parameters. All parameters must be printable,
/// and will actually be printed.
///
//// \brief an all-purpose exception class.
////////////////////////////////////////////////////////////////////////////
class GenericException{
public:
  //! O argument
  GenericException()
  {}
  //! 1 argument
  template<class A> GenericException(A x) 
    {
      cout<<endl<<endl<<x<<endl<<endl;
    }
  //! 2 arguments
  template<class A,class B> GenericException(A x,B y) 
    {
      cout<<endl<<endl<<x<<" "<<y<<endl<<endl;
    }
  //! 3 arguments
  template<class A,class B,class C> GenericException(A x,B y,C z) 
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<endl<<endl;
    }
  //! 4 arguments
  template<class A,class B,class C,class D> GenericException(A x,B y,C z,D a)
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<" "<<a<<endl<<endl;
    }
  //! 5 arguments
  template<class A,class B,class C,class D,class E> 
    GenericException(A x,B y,C z,D a,E b)
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<" "<<a<<" "<<b<<endl<<endl;
    }
  //! 6 arguments
  template<class A,class B,class C,class D,class E,class F> 
    GenericException(A x,B y,C z,D a,E b,F c)
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<" "<<a<<" "<<b<<
	" "<<c<<endl<<endl;
    }
  //! 7 arguments:
  template<class A,class B,class C,class D,class E,class F,class G> 
  GenericException(A x,B y,C z,D a,E b,F c,G d)
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<" "<<a<<" "<<b<<
	" "<<c<<" "<<d<<endl<<endl;
    }
  //! 8 arguments:
  template<class A,class B,class C,class D,class E,class F,class G,class H> 
  GenericException(A x,B y,C z,D a,E b,F c,G d,H e)
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<" "<<a<<" "<<b<<
	" "<<c<<" "<<d<<" "<<e<<endl<<endl;
    }
  //! 9 arguments:
  template<class A,class B,class C,class D,class E,class F,class G,class H,
	   class I> 
  GenericException(A x,B y,C z,D a,E b,F c,G d,H e,I f)
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<" "<<a<<" "<<b<<
	" "<<c<<" "<<d<<" "<<e<<" "<<f<<endl<<endl;
    }
  //! 10 arguments:
  template<class A,class B,class C,class D,class E,class F,class G,class H,
	   class I,class J> 
  GenericException(A x,B y,C z,D a,E b,F c,G d,H e,I f,J g)
    {
      cout<<endl<<endl<<x<<" "<<y<<" "<<z<<" "<<a<<" "<<b<<
	" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<g<<endl<<endl;
    }
};
#endif
