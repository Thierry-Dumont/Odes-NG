// define the "restrict" keyword for different compilers.
// define a ASSUME_IS_ALIGNED macro for all compilers.
// In AllocateDestroyVector.hpp, this is turned in ASSUME_ALIGNED, depending
// on compiler options.
//
// If the macro is not defined for a given compiler, a compilation error
// will occur.
#define BLOCK_FACTOR 8*sizeof(double)
#ifdef GCC
// gcc:
#define Restrict __restrict__
#define ASSUME_IS_ALIGNED(lvalueptr, align) lvalueptr = ( __typeof(lvalueptr))  __builtin_assume_aligned(lvalueptr, align)
#define ASSUME(l);
#define prag(l) GCC l
#endif
#ifdef ICC
// Intel compiler:
#define Restrict restrict
#define ASSUME_IS_ALIGNED(lvalueptr, align) __assume_aligned(lvalueptr, align)
#define ASSUME(L) __assume(L);
#define prag(l) l
#endif
#ifdef CLANG
// clang++
#define Restrict __restrict__
#define ASSUME_IS_ALIGNED(lvalueptr, align) 
#define ASSUME(l);
#endif
