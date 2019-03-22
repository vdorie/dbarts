#include "config.h"
#include <misc/simd.h>

#include <stdbool.h>
#include <misc/stddef.h>
#include <stdint.h>

#include <misc/types.h>
#include <misc/intrinsic.h>

/* static const char* const simdNames[] = {
  "none",
  "SSE",
  "SSE2",
  "SSE3",
  "SSSE3",
  "SSE4.1",
  "SSE4.2",
  "AVX",
  "AVX2",
  "AVX512F",
  "AVX512VL",
  "AVX512BW",
  "invalid"
}; */

#if !defined(__i386) && !defined(_X86_) && !defined(__x86_64__) && !defined(_M_AMD64) && !defined (_M_X64)
static misc_simd_instructionSet misc_simd_getMaxSIMDInstructionSet(void)
  return MISC_INST_C;
}
#else

#ifdef __GNUC__
#  include <cpuid.h>
#  ifndef __clang__
#    include <x86intrin.h>
#    include <immintrin.h>
#  endif
#elif defined(_MSC_VER)
#  include <intrin.h>
#elif defined(__INTEL_COMPILER)
#  include <immintrin.h>
#endif

static inline uint64_t xgetbv() {
#if (defined(_MSC_VER) && _MSC_FULL_VER >= 160040000) || \
    (defined(__INTEL_COMPILER) && __INTEL_COMPILER >= 1200) || \
    (defined(__XSAVE__) && defined(__GNUC__) && (__GNUC__ > 8 || (__GNUC__ == 8 && __GNUC_MINOR__ >= 2)))

  return (uint64_t) _xgetbv(0);

#elif defined(__GNUC__) || defined(__SUNPRO_C)
  
  uint32_t eax, edx;
  __asm__ (".byte 0x0f, 0x01, 0xd0" : "=a" (eax), "=d"(edx) : "c" (0));
  return ((uint64_t) edx << 32) | eax;

#else

  uint32_t veax, vedx;
  __asm {
    xor ecx, ecx    // for xcr0, set ecx to 0
    _emit 0x0f
    _emit 0x01
    _emit 0xd0
    mov veax, eax
    mov vedx, edx
  }
  return ((uint64_t) vedx << 32) | veax;
#endif
}

static inline int cpuid(int leafNumber, int* info)
{
#if defined(__GNUC__) || defined(__clang__)
  unsigned int* u_info =  (unsigned int*) info;
  return __get_cpuid(leafNumber, u_info, u_info + 1, u_info + 2, u_info + 3);
#elif defined (_MSC_VER) || defined (__INTEL_COMPILER)
  __cpuid(info, leafNumber);
#elif defined(__SUNPRO_C)
  __asm__ __volatile__(
    // check for 64 bit system
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
    "pushq %%rbx        \n\t" /* save %rbx */
#else
    "pushl %%ebx        \n\t" /* save %ebx */
#endif
    "cpuid              \n\t"
    "movl %%ebx ,%[ebx] \n\t" /* write the result into output var */
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
    "popq %%rbx         \n\t"
#else
    "popl %%ebx         \n\t"
#endif
    : "=a"(info[0]), [ebx] "=r"(info[1]), "=c"(info[2]), "=d"(info[3])
    : "a"(leafNumber));
#else
  __asm {
    mov eax, leafNumber
    xor ecx, ecx
    cpuid;
    mov esi, info
    mov [esi],      eax
    mov [esi + 4],  ebx
    mov [esi + 8],  ecx
    mov [esi + 12], edx
  }
#endif
  return info[0]; 
}

misc_simd_instructionSet misc_simd_getMaxSIMDInstructionSet(void)
{
  static misc_simd_instructionSet instructionSet = MISC_INST_INVALID;
  
  if (instructionSet != MISC_INST_INVALID) return instructionSet;
  
  instructionSet = MISC_INST_C;
  
  unsigned int maxBasicLeaf;
  int info[4] = { 0, 0, 0, 0 };
  
#if defined(__GNUC__)
  if ((maxBasicLeaf = __get_cpuid_max(0, NULL)) == 0) return instructionSet;
#else
  if (cpuid(0, info) == 0) return instructionSet;
  maxBasicLeaf = info[0];
#endif
  
  cpuid(1, info);
  
  int ecx = info[2], edx = info[3];
  
  if ((edx & (1 <<  0)) == 0) return instructionSet; // no floating point
  if ((edx & (1 << 23)) == 0) return instructionSet; // no MMX
  if ((edx & (1 << 15)) == 0) return instructionSet; // no conditional move
  if ((edx & (1 << 24)) == 0) return instructionSet; // no FXSAVE
  if ((edx & (1 << 25)) == 0) return instructionSet; // no SSE
  instructionSet = MISC_INST_SSE;
  
  if ((edx & (1 << 26)) == 0) return instructionSet; // no SSE2
  instructionSet = MISC_INST_SSE2;
  if ((ecx & (1 <<  0)) == 0) return instructionSet; // no SSE3
  instructionSet = MISC_INST_SSE3;
  if ((ecx & (1 <<  9)) == 0) return instructionSet; // no SSSE3
  instructionSet = MISC_INST_SSSE3;
  if ((ecx & (1 << 19)) == 0) return instructionSet; // no SSE4.1
  instructionSet = MISC_INST_SSE4_1;
  if ((ecx & (1 << 23)) == 0) return instructionSet; // no POPCNT
  if ((ecx & (1 << 20)) == 0) return instructionSet; // no SSE4.2
  instructionSet = MISC_INST_SSE4_2;
  if ((ecx & (1 << 27)) == 0) return instructionSet; // no OSXSAVE
  if ((xgetbv() & 6) != 6)    return instructionSet; // AVX not enabled in O.S.
  if ((ecx & (1 << 28)) == 0) return instructionSet; // no AVX
  instructionSet = MISC_INST_AVX;
  
  // check next leaf for feature flags
  if (maxBasicLeaf < 7) return instructionSet;
  
  cpuid(7, info);
  
  int ebx = info[1];  
  
  if ((ebx & (1 <<  5)) == 0) return instructionSet; // no AVX2
  instructionSet = MISC_INST_AVX2;
  if ((ebx & (1 << 16)) == 0) return instructionSet; // no AVX512
  
  if (maxBasicLeaf < 0x60) return instructionSet;
  
  cpuid(0xD, info);
  if ((info[0] & 0x60) != 0x60) return instructionSet; // no AVX512F
  instructionSet = MISC_INST_AVX512F; 
  
  if ((ebx & 0x80000000) == 0) return instructionSet; // no AVX512VL
  instructionSet = MISC_INST_AVX512VL; 
  
  if ((ebx & 0x40020000) != 0x40020000) return instructionSet; // no AVX512BW, AVX512DQ
  instructionSet = MISC_INST_AVX512BW; 
  
  return instructionSet;
}
#endif // !defined(__i386) && !defined(_X86_) && !defined(__x86_64__)

/* function pointers */
// partition
extern size_t (*misc_partitionRange)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t (*misc_partitionIndices)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
// linear algebra
extern void (*misc_addVectors)(const double* restrict x, misc_size_t length, double alpha, const double* restrict y, double* restrict z);
extern void (*misc_setVectorToConstant)(double* x, size_t length, double alpha);


/* implementing functions */
#ifdef HAVE_AVX2
extern size_t misc_partitionRange_avx2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_avx2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
#endif

#ifdef HAVE_AVX
extern void misc_addVectors_avx(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z);
extern void misc_setVectorToConstant_avx(double* x, size_t length, double alpha);
#endif

#ifdef HAVE_SSE4_1
extern size_t misc_partitionRange_sse4_1(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_sse4_1(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
#endif

#ifdef HAVE_SSE2
extern size_t misc_partitionRange_sse2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_sse2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern void misc_addVectors_sse2(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z);
extern void misc_setVectorToConstant_sse2(double* x, size_t length, double alpha);
#endif

// partition
extern size_t misc_partitionRange_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
// linearAlgebra
extern void misc_addVectors_c(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z);
extern void misc_setVectorToConstant_c(double* x, size_t length, double alpha);

void misc_simd_init(void) {
  misc_simd_instructionSet i = misc_simd_getMaxSIMDInstructionSet();
  
  misc_simd_setSIMDInstructionSet(i);
}

extern void misc_stat_setSIMDInstructionSet(misc_simd_instructionSet i);

void misc_simd_setSIMDInstructionSet(misc_simd_instructionSet i)
{
  // if the compiler supports SIMD instruction sets then the binary will be built
  // with code that supports them; check this by preprocessor defines
  
  // if the processor supports the SIMD instructions sets, set function pointers
  // to those entry points; check this by using cpuid
  
  // default to base C implementations
  
  if (i < MISC_INST_C || i >= MISC_INST_INVALID) return;
  
/* #ifdef HAVE_AVX2
#  define MAX_SIMD 8
#elif defined(HAVE_AVX)
#  define MAX_SIMD 7
#elif defined(HAVE_SSE4_2)
#  define MAX_SIMD 6
#elif defined(HAVE_SSE4_1)
#  define MAX_SIMD 5
#elif defined(HAVE_SSE3)
#  define MAX_SIMD 3
#elif defined(HAVE_SSE2)
#  define MAX_SIMD 2
#elif defined(HAVE_SSE)
#  define MAX_SIMD 1
#else
#  define MAX_SIMD 0
#endif */
  misc_simd_instructionSet i_max = misc_simd_getMaxSIMDInstructionSet();
  //ext_printf("SIMD instruction set %s (%d) requested, %s max compiled, %s max supported\n", simdNames[i], i,
  //           simdNames[MAX_SIMD], simdNames[i_max]);
  if (i > i_max) i = i_max;
  
#ifdef HAVE_AVX2
  if (i >= MISC_INST_AVX2) {
    misc_partitionRange = &misc_partitionRange_avx2;
    misc_partitionIndices = &misc_partitionIndices_avx2;
  } else
#endif
#ifdef HAVE_SSE4_1
  if (i >= MISC_INST_SSE4_1) {
    misc_partitionRange = &misc_partitionRange_sse4_1;
    misc_partitionIndices = &misc_partitionIndices_sse4_1;
  } else
#endif
#ifdef HAVE_SSE2
  if (i >= MISC_INST_SSE2) {
    misc_partitionRange = &misc_partitionRange_sse2;
    misc_partitionIndices = &misc_partitionIndices_sse2;
  } else
#endif
  {
    misc_partitionRange = &misc_partitionRange_c;
    misc_partitionIndices = &misc_partitionIndices_c;
  }
  
#ifdef HAVE_AVX
  if (i >= MISC_INST_AVX) {
    misc_addVectors = &misc_addVectors_avx;
    misc_setVectorToConstant = &misc_setVectorToConstant_avx;
  } else 
#endif
#ifdef HAVE_SSE2
  if (i >= MISC_INST_SSE2) {
    misc_addVectors = &misc_addVectors_sse2;
    misc_setVectorToConstant = &misc_setVectorToConstant_sse2;
  } else
#endif
  {
    misc_addVectors = &misc_addVectors_c;
    misc_setVectorToConstant = &misc_setVectorToConstant_c;
  }
  
  misc_stat_setSIMDInstructionSet(i);
}

