/* -*- C -*-
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2025  The R Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.

   It is part of the API and supports 'standalone Rmath'.
   Some entries possibly are not yet documented in 'Writing R Extensions'.


 * Reduced for dbarts from R 4.6.2's installed Rmath.h: only the standalone
 * branch and the declarations the copied files need are kept.
*/
#ifndef RMATH_H
#define RMATH_H

/* needed for cospi etc */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif

#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
# include <cmath>
#else
# include <math.h>
#endif

#ifdef NO_C_HEADERS
# warning "use of NO_C_HEADERS is defunct and will be ignored"
#endif

/*-- Mathlib as part of R --  define this for standalone : */
/* #undef MATHLIB_STANDALONE */

#define R_VERSION_STRING "4.6.2"

// Legacy defines -- C99 functions which R >= 3.5.0 requires
#ifndef HAVE_EXPM1
# define HAVE_EXPM1 1
#endif
#ifndef HAVE_HYPOT
# define HAVE_HYPOT 1
#endif
#ifndef HAVE_LOG1P
# define HAVE_LOG1P 1
#endif

#ifndef HAVE_WORKING_LOG1P
# define HAVE_WORKING_LOG1P 1
#endif

#if !defined(HAVE_WORKING_LOG1P)
/* remap to avoid problems with getting the right entry point */
extern "C" {
double  Rlog1p(double);
}
#define log1p Rlog1p
#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif


#ifdef MATHLIB_STANDALONE
# ifndef R_EXT_BOOLEAN_H_
/* "copy-paste" R_ext/Boolean.h if not already included: 
   This is standalone, so we do not worry about the size
*/
 #define R_EXT_BOOLEAN_H_
 #undef FALSE
 #undef TRUE
 typedef enum { FALSE = 0, TRUE } Rboolean;
 // stdbool.h is obsolescent but we could pro tem just use it.
 // But to future-proof copy from Boolean.h
 #if defined __STDC_VERSION__ && __STDC_VERSION__ > 202000L
 #elif defined __cplusplus
 #else
  # include <stdbool.h>
 #endif

# endif
#else
# error "this reduced copy of Rmath.h serves MATHLIB_STANDALONE only"
#endif


#define dnorm dnorm4
#define pnorm pnorm5
#define qnorm qnorm5

#ifdef  __cplusplus
extern "C" {
#endif
	/* Normal Distribution */

double	dnorm(double, double, double, int);
double	pnorm(double, double, double, int, int);
double	qnorm(double, double, double, int, int);
void	pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Gamma Distribution */

double	dgamma(double, double, double, int);
double	pgamma(double, double, double, int, int);
double	qgamma(double, double, double, int, int);

double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */

/* Compute the log of a sum or difference from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 * or  log (exp (logx) - exp (logy))
 *
 * without causing overflows or throwing away too much accuracy:
 */
double  logspace_add(double logx, double logy);
double  logspace_sub(double logx, double logy);
double  logspace_sum(const double *, int);

	/* Beta Distribution */

double	pbeta(double, double, double, int, int);

	/* Chi-squared Distribution */

double	qchisq(double, double, int, int);

	/* Student t Distribution */

double	dt(double, double, int);
double	pt(double, double, int, int);
double	qt(double, double, int, int);

	/* Poisson Distribution */

double	dpois_raw (double, double, int);
double	dpois(double, double, int);

	/* Gamma and Related Functions */
double	gammafn(double);
double	lgammafn(double);
double	lgammafn_sign(double, int*);

double	lbeta(double, double);

double	choose(double, double);
double	lchoose(double, double);

	/* General Support Functions */

double	fmax2(double, double);
double	fmin2(double, double);

/* More accurate cos(pi*x), sin(pi*x), tan(pi*x)

   These declarations might clash with system headers if someone had
   already included math.h with __STDC_WANT_IEC_60559_FUNCS_EXT__
   defined (and we try, above).
   We check for that via the value of __STDC_IEC_60559_FUNCS__
*/
#if !(defined(__STDC_IEC_60559_FUNCS__) && __STDC_IEC_60559_FUNCS__ >= 201506L)
double cospi(double);
double sinpi(double);
double tanpi(double);
#endif
double Rtanpi(double); /* our own in any case */

#ifdef  __cplusplus
}
#endif

#endif /* RMATH_H */
