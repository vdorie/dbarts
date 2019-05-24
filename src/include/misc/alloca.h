#ifndef MISC_ALLOCA_H
#define MISC_ALLOCA_H

#include <misc/stddef.h>

// Allocates on stack, if alloca exists if not, uses malloc. Thus, a free
// call is required but hopefully does nothing.
//
// define HAVE_ALLOCA_H to use

#ifdef __cplusplus
extern "C" {
#endif
  
#ifdef HAVE_ALLOCA_H
#  include <alloca.h>
#elif !defined alloca // try to pull something into scope
#  ifdef __GNUC__
#    define alloca __builtin_alloca
#  elif defined(__DECC)
#    define alloca __ALLOCA
#  elif defined(_MSC_VER)
#    include <malloc.h>
#    define alloca _alloca
#  else
#    ifdef __cplusplus
#      include <cstdlib>
#    else
#      include <stdlib.h>
#    endif
#  endif
#endif
  
#ifdef __cplusplus
}
#endif

#ifdef alloca
#  ifdef __cplusplus
#    define misc_stackAllocate(_N_, _T_) static_cast<_T_*>(alloca(static_cast<misc_size_t>(_N_) * sizeof(_T_)))
#  else
#    define misc_stackAllocate(_N_, _T_) (_T_*) alloca(((misc_size_t) (_N_)) * sizeof(_T_))
#  endif
#  define misc_stackFree(_P_)
#else
#  ifdef __cplusplus
#    define misc_stackAllocate(_N_, _T_) static_cast<_T_*>(::operator new(static_cast<misc_size_t>(_N_) * sizeof(_T_)))
#    define misc_stackFree(_P_) ::operator delete(_P_)
#  else
#    define misc_stackAllocate(_N_, _T_) (_T_*) malloc(((misc_size_t) (_N_)) * sizeof(_T_))
#    define misc_stackFree(_P_) free(_P_)
#  endif
#endif

#endif // MISC_ALLOCA_H

