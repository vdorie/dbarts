#include <misc/stats.h>

// lgamma stores the sign of Gamma(x) in the process-wide signgam, a data race
// when two threads call it; its reentrant form takes the sign by pointer and
// computes the same value, lgamma being that form over &signgam in every C
// library that has both. The reentrant form is in neither ISO C nor POSIX, so
// <math.h> declares it only under a feature macro that differs by C library
// (_REENTRANT on macOS, _DEFAULT_SOURCE on glibc), and it is declared here
// instead. MinGW-w64, whose libmingwex supplies the lgamma Rtools links in
// place of UCRT's, which has no reentrant form, names it __lgamma_r and
// declares it in no header.

#ifdef __MINGW32__

extern double __lgamma_r(double x, int* sign);

double misc_logGamma(double x)
{
  int sign;
  return __lgamma_r(x, &sign);
}

#else

extern double lgamma_r(double x, int* sign);

double misc_logGamma(double x)
{
  int sign;
  return lgamma_r(x, &sign);
}

#endif
