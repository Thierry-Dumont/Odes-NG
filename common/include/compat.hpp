#ifndef compat__h
#define compat__h
// We want to catch at compile time the erro:  Hessenberg=true and full=false.
//
// We do not want do depend of  something like boost::static_assert, and
// the compiler is supposed not to support static_assertions.
// This class will certainly diseapear when static_assert will be available
// on most compilers.

// Here, we use the trick defined at:
//   http://www.pixelbeat.org/programming/gcc/static_assert.html
// The compiler error message is not very clear, but it works.

// ---> start of macros defined at   
//http://www.pixelbeat.org/programming/gcc/static_assert.html:

//Note we need the 2 concats below because arguments to ##
// are not expanded, so we need to expand __LINE__ with one indirection
// before doing the actual concatenation. 
#define ASSERT_CONCAT_(a, b) a##b
#define ASSERT_CONCAT(a, b) ASSERT_CONCAT_(a, b)
#define ct_assert(e) enum { ASSERT_CONCAT(assert_line_, __LINE__) = 1/(!!(e)) }
// <--- end of macros.

template<bool full,bool Hessenberg> struct compat
{
  ct_assert(full || (!full && !Hessenberg));
};


#endif
