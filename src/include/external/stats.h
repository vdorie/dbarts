#ifndef EXTERNAL_STATS_H
#define EXTERNAL_STATS_H

#ifndef R_NO_REMAP_RMATH
#  define UNMAP_R_NO_REMAP_RMATH
#  define R_NO_REMAP
#  define R_NO_REMAP_RMATH
#endif

#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#define USE_FC_LEN_T
#endif

#include <R.h>
#include <Rmath.h> // used to pull in qchisq, et al

#undef USE_FC_LEN_T

#ifdef __cplusplus
extern "C" {
#endif

#define ext_quantileOfChiSquared(_P_, _NU_) Rf_qchisq((_P_), (_NU_), 1, 0)
#define ext_percentileOfChiSquared(_Q_, _NU_) Rf_pchisq((_Q_), (_NU_), 1, 0)
  
#define ext_densityOfNormal(_X_, _MU_, _SIGMA_) Rf_dnorm4((_X_), (_MU_), (_SIGMA_), 0)
#define ext_cumulativeProbabilityOfNormal(_Q_, _MU_, _SIGMA_) Rf_pnorm5((_Q_), (_MU_), (_SIGMA_), 1, 0)
#define ext_quantileOfNormal(_P_, _MU_, _SIGMA_) Rf_qnorm5((_P_), (_MU_), (_SIGMA_), 1, 0)

extern double Rf_pnorm5(double, double, double, int, int);
extern double Rf_pchisq(double, double, int, int);
extern double Rf_qchisq(double, double, int, int);

#ifdef __cplusplus
}
#endif

#ifdef UNMAP_R_NO_REMAP_RMATH
#  undef R_NO_REMAP_RMATH
#  undef R_NO_REMAP
#  undef UNMAP_R_NO_REMAP_RMATH
#endif

#endif // EXTERNAL_STATS_H

