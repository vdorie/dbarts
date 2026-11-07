#ifndef EXTERNAL_STATS_H
#define EXTERNAL_STATS_H

#ifndef DBARTS_HOST_STANDALONE

#ifndef R_NO_REMAP_RMATH
#  define UNMAP_R_NO_REMAP_RMATH
#  ifndef R_NO_REMAP
#    define UNMAP_R_NO_REMAP
#    define R_NO_REMAP
#  endif
#  define R_NO_REMAP_RMATH
#endif

#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#define USE_FC_LEN_T
#endif

#include <Rmath.h> // used to pull in qchisq, et al

#undef USE_FC_LEN_T

#ifdef UNMAP_R_NO_REMAP_RMATH
#  undef R_NO_REMAP_RMATH
#  ifdef UNMAP_R_NO_REMAP
#    undef R_NO_REMAP
#    undef UNMAP_R_NO_REMAP
#  endif
#  undef UNMAP_R_NO_REMAP_RMATH
#endif

#endif // DBARTS_HOST_STANDALONE

#include <external/nmath.h> // the seven ext_ names and their implementation

#define ext_quantileOfChiSquared(_P_, _NU_) ext_qchisq((_P_), (_NU_), 1, 0)

#define ext_densityOfNormal(_X_, _MU_, _SIGMA_) ext_dnorm((_X_), (_MU_), (_SIGMA_), 0)
#define ext_cumulativeProbabilityOfNormal(_Q_, _MU_, _SIGMA_) ext_pnorm((_Q_), (_MU_), (_SIGMA_), 1, 0)
#define ext_quantileOfNormal(_P_, _MU_, _SIGMA_) ext_qnorm((_P_), (_MU_), (_SIGMA_), 1, 0)
#define ext_quantileOfStudentT(_P_, _NU_) ext_qt((_P_), (_NU_), 1, 0)

#endif // EXTERNAL_STATS_H
