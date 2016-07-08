#ifndef AllocateDestroyVector__h
#define AllocateDestroyVector__h
//////////////////////////////////////////////////////////////////////
/// Define different methods to allocate arrays.
//////////////////////////////////////////////////////////////////////
#include "MacrosForCompilers.hpp"
#ifdef ALIGN_BLOCKED
#ifdef ICC
#include <malloc.h>
#else 
#include "mm_malloc.h"
#endif
// Blocked alignment allocation.
double *allocDoubleArray(int size)
{return (double *) _mm_malloc(size*sizeof(double),BLOCK_FACTOR);}
void destroyDoubleArray(double *x){_mm_free(x);}
#define ASSUME_ALIGNED(lvalueptr) ASSUME_IS_ALIGNED(lvalueptr,BLOCK_FACTOR)
#else
// Classical {new,delete} allocation. We cannot assume aligment.
double *allocDoubleArray(int size){return new double[size];}
void destroyDoubleArray(double *x){delete[] x;}
#define ASSUME_ALIGNED(lvalueptr)
#endif
#endif
