#ifndef EXTERNAL_ALLOCA_H
#define EXTERNAL_ALLOCA_H

#include "stddef.h"

// Allocates on stack, if alloca exists if not, uses malloc. Thus, a free
// call is required but hopefully does nothing.
//
// define HAVE_ALLOCA_H to use

#ifdef __cplusplus
extern "C" {
#endif
  
#ifdef HAVE_ALLOCA_H
#  include <alloca.h>
#elif !defined alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# elif defined(__DECC)
#  define alloca __ALLOCA
# elif defined(_MSC_VER)
#  include <malloc.h>
#  define alloca _alloca
# elif defined(__sun)
#  include <alloca.h>
# else
#  ifdef __cplusplus
#   include <cstdlib>
#  else
#   include <stdlib.h>
#  endif
# endif
#endif
  
#ifdef __cplusplus
}
#endif

#ifdef alloca
# define ext_stackAllocate(_N_, _T_) (_T_ *) alloca(((ext_size_t) (_N_)) * sizeof(_T_))
# define ext_stackFree(_P_)
#else
# define ext_stackAllocate(_N_, _T_) (_T_ *) malloc(((ext_size_t) (_N_)) * sizeof(_T_))
# define ext_stackFree(_P_) free(_P_)
#endif

#endif // EXTERNAL_ALLOCA_H
