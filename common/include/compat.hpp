#ifndef compat__h
#define compat__h
#include <type_traits>
//
// We want to catch at compile time the error:  Hessenberg=true and full=false.
//
template<bool full,bool Hessenberg> struct compat
{
  static_assert(full || (!full && !Hessenberg),
	  "\n\n Hessenberg option can only be used with full Jacobian.\n\n");
};


#endif
