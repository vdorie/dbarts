// file includes global variables

#ifndef R_INTERFACE_HPP
#define R_INTERFACE_HPP

#ifdef THREAD_SAFE_UNLOAD
#  include <pthread.h>
#endif

#include <set>

// struct SEXPREC;
// typedef struct SEXPREC* SEXP;
#include <external/Rinternals.h> // SEXP

typedef bool(*ExternalPointerComparator)(const SEXP&lhs, const SEXP& rhs);

typedef std::set<SEXP, ExternalPointerComparator> PointerSet;
extern PointerSet* activeFits;

#ifdef THREAD_SAFE_UNLOAD
extern pthread_mutex_t fitMutex;
#endif

#endif

