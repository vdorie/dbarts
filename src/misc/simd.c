#include "config.h"
#include <misc/simd.h>

#include <misc/stddef.h>
#include <misc/types.h>
#include <misc/intrinsic.h>

#define SIMD_SSE      0x1
#define SIMD_SSE2     0x2
#define SIMD_SSE3     0x4
#define SIMD_SSE4_1   0x8
#define SIMD_SSE4_2  0x10
#define SIMD_AVX     0x20
#define SIMD_AVX2    0x40
#define SIMD_AVX512F 0x80
 
static unsigned int getx86SIMDSets(void) {
#ifdef __SUNPRO_C
  return 0;
#else
  unsigned int eax, ebx, ecx, edx, flag = 0;
#  ifdef _MSC_VER
  int cpuid[4];
  __cpuid(cpuid, 1);
  eax = cpuid[0], ebx = cpuid[1], ecx = cpuid[2], edx = cpuid[3];
#  elif __GNUC__
  __asm__("cpuid" : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "a" (1));
#  else
  asm volatile("cpuid" : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "a" (1));
#  endif
  
  if ((edx >> 25) & 1) flag |= SIMD_SSE;
  if ((edx >> 26) & 1) flag |= SIMD_SSE2;
  if ((ecx >>  0) & 1) flag |= SIMD_SSE3;
  if ((ecx >> 19) & 1) flag |= SIMD_SSE4_1;
  if ((ecx >> 20) & 1) flag |= SIMD_SSE4_2;
  if ((ecx >> 28) & 1) flag |= SIMD_AVX;
  if ((ebx >>  5) & 1) flag |= SIMD_AVX2;
  if ((ebx >> 16) & 1) flag |= SIMD_AVX512F;
  
  return flag;
#endif
}

misc_simd_instructionLevel misc_simd_getMaxSIMDInstructionSet() {
  unsigned int flag = getx86SIMDSets();
  
  return flag == 0 ? MISC_INST_C : (misc_simd_instructionLevel) (8 * sizeof(unsigned int) - countLeadingZeros(flag));
}

// function pointers
extern size_t (*misc_partitionRange)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t (*misc_partitionIndices)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern void (*misc_addVectors)(const double* restrict x, misc_size_t length, double alpha, const double* restrict y, double* restrict z);
extern void (*misc_setVectorToConstant)(double* x, size_t length, double alpha);

// implementing functions
extern size_t misc_partitionRange_avx2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionRange_sse4_1(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionRange_sse2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionRange_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_avx2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_sse4_1(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_sse2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern void misc_addVectors_c(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z);
extern void misc_addVectors_sse2(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z);
extern void misc_addVectors_avx(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z);
extern void misc_setVectorToConstant_c(double* x, size_t length, double alpha);
extern void misc_setVectorToConstant_sse2(double* x, size_t length, double alpha);
extern void misc_setVectorToConstant_avx(double* x, size_t length, double alpha);

extern void misc_initStats(misc_simd_instructionLevel i);
extern void misc_initART(misc_simd_instructionLevel i);

void misc_simd_init(void) {
  int flag = getx86SIMDSets();
  
  // if the compiler supports SIMD instruction sets then the binary will be built
  // with code that supports them; check this by preprocessor defines
  
  // if the processor supports the SIMD instructions sets, set function pointers
  // to those entry points; check this by using cpuid
  
  // default to base C implementations
  
#ifdef __AVX2__
  if (flag & SIMD_AVX2) {
    misc_partitionRange = &misc_partitionRange_avx2;
    misc_partitionIndices = &misc_partitionIndices_avx2;
    misc_addVectors = &misc_addVectors_avx;
    misc_setVectorToConstant = &misc_setVectorToConstant_avx;
  } else 
#endif
#ifdef __AVX__
  if (flag & SIMD_AVX) {
    misc_partitionRange = &misc_partitionRange_sse4_1;
    misc_partitionIndices = &misc_partitionIndices_sse4_1;
    misc_addVectors = &misc_addVectors_avx;
    misc_setVectorToConstant = &misc_setVectorToConstant_avx;
  } else 
#endif
#ifdef __SSE4_1__
  if (flag & SIMD_SSE4_1) {
    misc_partitionRange = &misc_partitionRange_sse4_1;
    misc_partitionIndices = &misc_partitionIndices_sse4_1;
    misc_addVectors = &misc_addVectors_sse2;
    misc_setVectorToConstant = &misc_setVectorToConstant_sse2;
  } else 
#endif
#ifdef __SSE2__
  if (flag & SIMD_SSE2) {
    misc_partitionRange = &misc_partitionRange_sse2;
    misc_partitionIndices = &misc_partitionIndices_sse2;
    misc_addVectors = &misc_addVectors_sse2;
    misc_setVectorToConstant = &misc_setVectorToConstant_sse2;
  } else
#endif
  {
    misc_partitionRange = &misc_partitionRange_c;
    misc_partitionIndices = &misc_partitionIndices_c;
    misc_addVectors = &misc_addVectors_c;
    misc_setVectorToConstant = &misc_setVectorToConstant_c;
  }
  
  misc_simd_instructionLevel i = flag == 0 ? MISC_INST_C : (misc_simd_instructionLevel) (8 * sizeof(unsigned int) - countLeadingZeros(flag));
  
  misc_initART(i);
  misc_initStats(i);
}

// #include <external/io.h>

void misc_simd_setSIMDInstructionSet(misc_simd_instructionLevel i)
{
  if (i < MISC_INST_C || i > MISC_INST_AVX512F) return;
  
#ifdef __AVX2__
  if (i >= MISC_INST_AVX2) {
    // ext_printf("setting SIMD level to AVX2\n");
    misc_partitionRange = &misc_partitionRange_sse4_1;
    misc_partitionIndices = &misc_partitionIndices_sse4_1;
    misc_addVectors = &misc_addVectors_avx;
    misc_setVectorToConstant = &misc_setVectorToConstant_avx;
  } else
#endif
#ifdef __AVX__
  if (i >= MISC_INST_AVX) {
    // ext_printf("setting SIMD level to AVX\n");
    misc_partitionRange = &misc_partitionRange_sse4_1;
    misc_partitionIndices = &misc_partitionIndices_sse4_1;
    misc_addVectors = &misc_addVectors_avx;
    misc_setVectorToConstant = &misc_setVectorToConstant_avx;
  } else
#endif
#ifdef __SSE4_1__
  if (i >= MISC_INST_SSE4_1) {
    // ext_printf("setting SIMD level to SSE4_1\n");
    misc_partitionRange = &misc_partitionRange_sse4_1;
    misc_partitionIndices = &misc_partitionIndices_sse4_1;
    misc_addVectors = &misc_addVectors_sse2;
    misc_setVectorToConstant = &misc_setVectorToConstant_sse2;
  } else
#endif
#ifdef __SSE2__
  if (i >= MISC_INST_SSE2) {
    // ext_printf("setting SIMD level to SSE2\n");
    misc_partitionRange = &misc_partitionRange_sse2;
    misc_partitionIndices = &misc_partitionIndices_sse2;
    misc_addVectors = &misc_addVectors_sse2;
    misc_setVectorToConstant = &misc_setVectorToConstant_sse2;
  } else
#endif
  {
    // ext_printf("setting SIMD level to C\n");
    misc_partitionRange = &misc_partitionRange_c;
    misc_partitionIndices = &misc_partitionIndices_c;
    misc_addVectors = &misc_addVectors_c;
    misc_setVectorToConstant = &misc_setVectorToConstant_c;
  }
  
  misc_initART(i);
  misc_initStats(i);
}


