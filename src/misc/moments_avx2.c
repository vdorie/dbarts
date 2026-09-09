#include "config.h"

#include <misc/stddef.h>
#include <misc/stats.h>
#include <misc/intrinsic.h>

// Four-wide weighted suffstat kernels. Built only where the Makefile has an
// AVX2 flag to give it, and reached only through moments.c's function
// pointers, which misc_stat_setSIMDInstructionSet points here when the
// running CPU reports AVX2. The reference build never compiles this file's
// bodies: it has no vector suffstat kernel at all.
#if defined(COMPILER_SUPPORTS_AVX2) && !defined(DBARTS_REFERENCE_BUILD)

// One 256-bit register holds all four banks, lane b being bank b, which is
// EXACTLY the split the two-lane kernels in moments.c build out of two
// registers: the length % 4 prologue accumulates into bank 0, element i then
// goes to bank (i - length % 4) mod 4, and the combine is ((b0 + b1) + b2) + b3
// strictly left to right. Every sum below is therefore the same bytes this
// package computes on SSE2, on NEON and with no vector unit at all, which is
// what lets a runtime CPU check select it without moving a draw.
//
// Each product is named before it is accumulated: the sums must round w * x
// before adding it, so no fused multiply-add may swallow that rounding.

/// Left-to-right combine of the four banks a 256-bit accumulator holds.
static inline double combineBanks(__m256d banks)
{
  double bank[4];
  _mm256_storeu_pd(bank, banks);
  return ((bank[0] + bank[1]) + bank[2]) + bank[3];
}

void misc_computeWeightedSufficientStatistics_avx2(const double* restrict x, size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX)
{
  size_t i = 0, prologue = length % 4;
  double sw0 = 0.0, swx0 = 0.0;
  for ( ; i < prologue; ++i) { double wi = w[i], p = wi * x[i]; sw0 += wi; swx0 += p; }

  // bank 0 carries the prologue before any main-loop term reaches it
  __m256d sw = _mm256_set_pd(0.0, 0.0, 0.0, sw0);
  __m256d sx = _mm256_set_pd(0.0, 0.0, 0.0, swx0);
  for ( ; i < length; i += 4) {
    __m256d wv = _mm256_loadu_pd(w + i);
    __m256d p = _mm256_mul_pd(wv, _mm256_loadu_pd(x + i));
    sw = _mm256_add_pd(sw, wv);
    sx = _mm256_add_pd(sx, p);
  }

  *sumW = combineBanks(sw);
  *sumWX = combineBanks(sx);
}

void misc_computeIndexedWeightedSufficientStatistics_avx2(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX)
{
  size_t i = 0, prologue = length % 4;
  double sw0 = 0.0, swx0 = 0.0;
  for ( ; i < prologue; ++i) {
    size_t j = indices[i];
    double wi = w[j], p = wi * x[j];
    sw0 += wi; swx0 += p;
  }

  __m256d sw = _mm256_set_pd(0.0, 0.0, 0.0, sw0);
  __m256d sx = _mm256_set_pd(0.0, 0.0, 0.0, swx0);
  for ( ; i < length; i += 4) {
    // packed from four scalar loads on purpose: a hardware gather measured
    // slower than scalar loads on this kernel
    size_t j0 = indices[i], j1 = indices[i + 1], j2 = indices[i + 2], j3 = indices[i + 3];
    __m256d wv = _mm256_set_pd(w[j3], w[j2], w[j1], w[j0]);
    __m256d p = _mm256_mul_pd(wv, _mm256_set_pd(x[j3], x[j2], x[j1], x[j0]));
    sw = _mm256_add_pd(sw, wv);
    sx = _mm256_add_pd(sx, p);
  }

  *sumW = combineBanks(sw);
  *sumWX = combineBanks(sx);
}

#else

// ISO C forbids an empty translation unit, and every non-x86 build and every
// reference build compiles this file to exactly that.
typedef int momentsAvx2TranslationUnitNotEmpty;

#endif
