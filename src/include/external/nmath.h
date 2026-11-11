#ifndef EXTERNAL_NMATH_H
#define EXTERNAL_NMATH_H

// The six math routines the engine and the support library call, each with
// every flag spelled. Which implementation stands behind them is fixed when
// the code is compiled: R's own math library in the R package, so every draw
// is what R computes, and a copy of R's routines, renamed with a dbarts_nmath_
// prefix, when DBARTS_HOST_STANDALONE is defined. Neither branch includes
// Rmath.h.

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DBARTS_HOST_STANDALONE

extern double dbarts_nmath_dnorm4(double, double, double, int);
extern double dbarts_nmath_pnorm5(double, double, double, int, int);
extern double dbarts_nmath_qnorm5(double, double, double, int, int);
extern double dbarts_nmath_qt(double, double, int, int);
extern double dbarts_nmath_qchisq(double, double, int, int);
extern double dbarts_nmath_lchoose(double, double);

#  define ext_dnorm(_X_, _MU_, _SD_, _LOG_P_) dbarts_nmath_dnorm4((_X_), (_MU_), (_SD_), (_LOG_P_))
#  define ext_pnorm(_X_, _MU_, _SD_, _LOWER_TAIL_, _LOG_P_) \
     dbarts_nmath_pnorm5((_X_), (_MU_), (_SD_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_qnorm(_P_, _MU_, _SD_, _LOWER_TAIL_, _LOG_P_) \
     dbarts_nmath_qnorm5((_P_), (_MU_), (_SD_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_qt(_P_, _DF_, _LOWER_TAIL_, _LOG_P_) dbarts_nmath_qt((_P_), (_DF_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_qchisq(_P_, _DF_, _LOWER_TAIL_, _LOG_P_) \
     dbarts_nmath_qchisq((_P_), (_DF_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_lchoose(_N_, _K_) dbarts_nmath_lchoose((_N_), (_K_))

#else

extern double Rf_dnorm4(double, double, double, int);
extern double Rf_pnorm5(double, double, double, int, int);
extern double Rf_qnorm5(double, double, double, int, int);
extern double Rf_qt(double, double, int, int);
extern double Rf_qchisq(double, double, int, int);
extern double Rf_lchoose(double, double);

#  define ext_dnorm(_X_, _MU_, _SD_, _LOG_P_) Rf_dnorm4((_X_), (_MU_), (_SD_), (_LOG_P_))
#  define ext_pnorm(_X_, _MU_, _SD_, _LOWER_TAIL_, _LOG_P_) \
     Rf_pnorm5((_X_), (_MU_), (_SD_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_qnorm(_P_, _MU_, _SD_, _LOWER_TAIL_, _LOG_P_) \
     Rf_qnorm5((_P_), (_MU_), (_SD_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_qt(_P_, _DF_, _LOWER_TAIL_, _LOG_P_) Rf_qt((_P_), (_DF_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_qchisq(_P_, _DF_, _LOWER_TAIL_, _LOG_P_) Rf_qchisq((_P_), (_DF_), (_LOWER_TAIL_), (_LOG_P_))
#  define ext_lchoose(_N_, _K_) Rf_lchoose((_N_), (_K_))

#endif

#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_NMATH_H
