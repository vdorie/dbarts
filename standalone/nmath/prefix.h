/* Force-included ahead of every file of the copy (cc -include prefix.h).
 *
 * The first block stands in for R's config.h: the macros the copied files
 * read, set as R's configure sets them when it builds R itself, so the copy
 * compiles to R's own code. The second renames every external symbol the copy
 * defines, so a program that also links R or a standalone libRmath sees no
 * clash; the recipe that builds the copy fails if any defined symbol escapes
 * the prefix. */
#ifndef DBARTS_NMATH_PREFIX_H
#define DBARTS_NMATH_PREFIX_H

#define HAVE_LONG_DOUBLE 1
#define HAVE_NEARBYINT 1
#define HAVE_WORKING_ISFINITE 1
/* Apple's libm has __cospi and friends and no C23 cospi; elsewhere R's
   configure finds neither on the C libraries R supports today, and the copy
   uses its own. */
#if defined(__APPLE__)
#  define HAVE___COSPI 1
#  define HAVE___SINPI 1
#  define HAVE___TANPI 1
#endif

#define cospi                dbarts_nmath_cospi
#define d1mach_              dbarts_nmath_d1mach_
#define dgamma               dbarts_nmath_dgamma
#define dnorm4               dbarts_nmath_dnorm4
#define dpois                dbarts_nmath_dpois
#define dpois_raw            dbarts_nmath_dpois_raw
#define dt                   dbarts_nmath_dt
#define fmax2                dbarts_nmath_fmax2
#define fmin2                dbarts_nmath_fmin2
#define gammafn              dbarts_nmath_gammafn
#define i1mach_              dbarts_nmath_i1mach_
#define lbeta                dbarts_nmath_lbeta
#define lgamma1p             dbarts_nmath_lgamma1p
#define lgammafn             dbarts_nmath_lgammafn
#define lgammafn_sign        dbarts_nmath_lgammafn_sign
#define log1pmx              dbarts_nmath_log1pmx
#define logspace_add         dbarts_nmath_logspace_add
#define logspace_sub         dbarts_nmath_logspace_sub
#define logspace_sum         dbarts_nmath_logspace_sum
#define pbeta                dbarts_nmath_pbeta
#define pbeta_raw            dbarts_nmath_pbeta_raw
#define pgamma               dbarts_nmath_pgamma
#define pnorm_both           dbarts_nmath_pnorm_both
#define pnorm5               dbarts_nmath_pnorm5
#define pt                   dbarts_nmath_pt
#define qchisq               dbarts_nmath_qchisq
#define qchisq_appr          dbarts_nmath_qchisq_appr
#define qgamma               dbarts_nmath_qgamma
#define qnorm5               dbarts_nmath_qnorm5
#define qt                   dbarts_nmath_qt
#define Rf_bd0               dbarts_nmath_Rf_bd0
#define Rf_bratio            dbarts_nmath_Rf_bratio
#define Rf_chebyshev_eval    dbarts_nmath_Rf_chebyshev_eval
#define Rf_chebyshev_init    dbarts_nmath_Rf_chebyshev_init
#define Rf_d1mach            dbarts_nmath_Rf_d1mach
#define Rf_ebd0              dbarts_nmath_Rf_ebd0
#define Rf_i1mach            dbarts_nmath_Rf_i1mach
#define Rf_lgammacor         dbarts_nmath_Rf_lgammacor
#define Rf_pgamma_raw        dbarts_nmath_Rf_pgamma_raw
#define Rf_stirlerr          dbarts_nmath_Rf_stirlerr
#define Rtanpi               dbarts_nmath_Rtanpi
#define sinpi                dbarts_nmath_sinpi
#define tanpi                dbarts_nmath_tanpi

#endif /* DBARTS_NMATH_PREFIX_H */
