/* This is taken from bits and pieces all over the internet. The strongest inspiriation is
 * Agner Fog's C++ vector class library (https://www.agner.org/optimize/#vectorclass)
 * which is vailable under the GPL. The cross-platform cpuid is from Wikipedia 
 * (https://en.wikipedia.org/wiki/CPUID#CPUID_usage_from_high-level_languages) and is 
 * available under a creative-commons license. I've mostly tried to use built-ins
 * when possible, and probably done more testing on Solaris-x86 than most.
 */
#include "config.h"
#include <misc/simd.h>

#include <stdbool.h>
#include <misc/stddef.h>
#include <stdint.h>

#include <misc/types.h>
#include <misc/intrinsic.h>

int misc_simd_alignment = 0;

/* static const char* const simdNames[] = { // for x86
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

/* static const char* const simdNames[] = { // for arm
  "none",
  "neon",
  "sve",
  "sve2",
  "invalid"
}; */

// if not on any x86 descendent, use pure C no matter what
#if !defined(__i386) && !defined(_X86_) && !defined(__x86_64__) && !defined(_M_AMD64) && !defined (_M_X64)

#  if defined(__arm__) || defined(__aarch64__) || defined(_ARM) || defined(_M_ARM)
misc_simd_instructionSet misc_simd_getMaxSIMDInstructionSet(void) {
  // NOTE THAT THIS IS OBVIOUSLY NOT IMPLEMENTED YET. OS specific calls
  // are likely required, and it appears that SVE can be supported
  // with different vector lengths.
  return MISC_INST_NEON;
}
#  else
misc_simd_instructionSet misc_simd_getMaxSIMDInstructionSet(void) {
  return MISC_INST_C;
}
#  endif
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
    // check for 64 bit system, use rbx there and ebx on 32 bit
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
  if ((xgetbv() & 6) != 6)    return instructionSet; // AVX not enabled in OS
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
extern void (*misc_addVectorsInPlace)(const double* restrict x, misc_size_t length, double* restrict y);
extern void (*misc_subtractVectorsInPlace)(const double* restrict x, misc_size_t length, double* restrict y);
extern void (*misc_addVectorsInPlaceWithMultiplier)(const double* restrict x, misc_size_t length, double alpha, double* restrict y);

extern void (*misc_addAlignedVectorsInPlace)(const double* restrict x, misc_size_t length, double* restrict y);
extern void (*misc_subtractAlignedVectorsInPlace)(const double* restrict x, misc_size_t length, double* restrict y);

extern void (*misc_addScalarToVectorInPlace)(double* x, size_t length, double alpha);
extern void (*misc_setVectorToConstant)(double* x, misc_size_t length, double alpha);
extern void (*misc_transposeMatrix)(const double* restrict x, size_t numRows, size_t numCols, double* restrict y);


/* implementing functions */
#ifdef COMPILER_SUPPORTS_AVX2
extern size_t misc_partitionRange_avx2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_avx2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
#endif

#ifdef COMPILER_SUPPORTS_AVX
extern void misc_addVectorsInPlace_avx(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_subtractVectorsInPlace_avx(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_addVectorsInPlaceWithMultiplier_avx(const double* restrict x, misc_size_t length, double alpha, double* restrict y);

// extern void misc_addAlignedVectorsInPlace_avx(const double* restrict x, misc_size_t length, double* restrict y);
// extern void misc_subtractAlignedVectorsInPlace_avx(const double* restrict x, misc_size_t length, double* restrict y);
// extern void misc_addAlignedVectorsInPlaceWithMultiplier_avx(const double* restrict x, misc_size_t length, double alpha, double* restrict y);

extern void misc_addScalarToVectorInPlace_avx(double* x, misc_size_t length, double alpha);
extern void misc_setVectorToConstant_avx(double* x, misc_size_t length, double alpha);
extern void misc_transposeMatrix_avx(const double* restrict x, misc_size_t numRows, size_t numCols, double* restrict y);
#endif

#ifdef COMPILER_SUPPORTS_SSE4_1
extern size_t misc_partitionRange_sse4_1(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_sse4_1(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
#endif

#ifdef COMPILER_SUPPORTS_SSE2
extern size_t misc_partitionRange_sse2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_sse2(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);

extern void misc_addVectorsInPlace_sse2(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_subtractVectorsInPlace_sse2(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_addVectorsInPlaceWithMultiplier_sse2(const double* restrict x, misc_size_t length, double alpha, double* restrict y);

// extern void misc_addAlignedVectorsInPlace_sse2(const double* restrict x, misc_size_t length, double* restrict y);
// extern void misc_subtractAlignedVectorsInPlace_sse2(const double* restrict x, misc_size_t length, double* restrict y);

extern void misc_addScalarToVectorInPlace_sse2(double* x, misc_size_t length, double alpha);
extern void misc_setVectorToConstant_sse2(double* x, misc_size_t length, double alpha);
extern void misc_transposeMatrix_sse2(const double* restrict x, misc_size_t numRows, misc_size_t numCols, double* restrict y);
#endif

#ifdef COMPILER_SUPPORTS_NEON
extern size_t misc_partitionRange_neon(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_neon(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);

extern void misc_addVectorsInPlace_neon(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_subtractVectorsInPlace_neon(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_addVectorsInPlaceWithMultiplier_neon(const double* restrict x, misc_size_t length, double alpha, double* restrict y);

extern void misc_addAlignedVectorsInPlace_neon(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_subtractAlignedVectorsInPlace_neon(const double* restrict x, misc_size_t length, double* restrict y);

extern void misc_addScalarToVectorInPlace_neon(double* x, misc_size_t length, double alpha);
extern void misc_setVectorToConstant_neon(double* x, misc_size_t length, double alpha);
extern void misc_transposeMatrix_neon(const double* restrict x, misc_size_t numRows, misc_size_t numCols, double* restrict y);
#endif

// partition
extern size_t misc_partitionRange_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern size_t misc_partitionIndices_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
// linearAlgebra
extern void misc_addVectorsInPlace_c(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_subtractVectorsInPlace_c(const double* restrict x, misc_size_t length, double* restrict y);
extern void misc_addVectorsInPlaceWithMultiplier_c(const double* restrict x, misc_size_t length, double alpha, double* restrict y);

extern void misc_addScalarToVectorInPlace_c(double* x, misc_size_t length, double alpha);
extern void misc_setVectorToConstant_c(double* x, misc_size_t length, double alpha);
extern void misc_transposeMatrix_c(const double* restrict x, misc_size_t numRows, misc_size_t numCols, double* restrict y);

void misc_simd_init(void) {
  misc_simd_instructionSet i = misc_simd_getMaxSIMDInstructionSet();
  
  misc_simd_setSIMDInstructionSet(i);
}

extern void misc_stat_setSIMDInstructionSet(misc_simd_instructionSet i);

#include <external/io.h>

void misc_simd_setSIMDInstructionSet(misc_simd_instructionSet i)
{
  // if the compiler supports SIMD instruction sets then the binary will be built
  // with code that supports them; check this by preprocessor defines
  
  // if the processor supports the SIMD instructions sets, set function pointers
  // to those entry points; check this by using cpuid
  
  // default to base C implementations
  
  if (i < MISC_INST_C || i >= MISC_INST_INVALID) return;
  
  misc_simd_instructionSet i_max = misc_simd_getMaxSIMDInstructionSet();
  if (i > i_max) i = i_max;
  
  // Integer
#ifdef COMPILER_SUPPORTS_AVX2
  if (i >= MISC_INST_AVX2) {
    misc_partitionRange = &misc_partitionRange_avx2;
    misc_partitionIndices = &misc_partitionIndices_avx2;
  } else
#endif
#ifdef COMPILER_SUPPORTS_SSE4_1
  if (i >= MISC_INST_SSE4_1) {
    misc_partitionRange = &misc_partitionRange_sse4_1;
    misc_partitionIndices = &misc_partitionIndices_sse4_1;
  } else
#endif
#ifdef COMPILER_SUPPORTS_SSE2
  if (i >= MISC_INST_SSE2) {
    misc_partitionRange = &misc_partitionRange_sse2;
    misc_partitionIndices = &misc_partitionIndices_sse2;
  } else
#endif
#ifdef COMPILER_SUPPORTS_NEON
  if (i >= MISC_INST_NEON) {
    misc_partitionRange = &misc_partitionRange_neon;
    misc_partitionIndices = &misc_partitionIndices_neon;
  } else
#endif
  {
    misc_partitionRange = &misc_partitionRange_c;
    misc_partitionIndices = &misc_partitionIndices_c;
  }
  
  // Float
#ifdef COMPILER_SUPPORTS_AVX
  if (i >= MISC_INST_AVX) {
    misc_simd_alignment = 0;
    misc_addAlignedVectorsInPlace = &misc_addVectorsInPlace_avx;
    misc_subtractAlignedVectorsInPlace = &misc_subtractVectorsInPlace_avx;
    
    misc_addVectorsInPlace = &misc_addVectorsInPlace_avx;
    misc_subtractVectorsInPlace = &misc_subtractVectorsInPlace_avx;
    misc_addVectorsInPlaceWithMultiplier = &misc_addVectorsInPlaceWithMultiplier_avx;
   
    misc_addScalarToVectorInPlace = &misc_addScalarToVectorInPlace_avx;
    misc_setVectorToConstant = &misc_setVectorToConstant_avx;
    misc_transposeMatrix = &misc_transposeMatrix_avx;
  } else 
#endif
#ifdef COMPILER_SUPPORTS_SSE2
  if (i >= MISC_INST_SSE2) {
    misc_simd_alignment = 0;
    misc_addAlignedVectorsInPlace = &misc_addVectorsInPlace_sse2;
    misc_subtractAlignedVectorsInPlace = &misc_subtractVectorsInPlace_sse2;
    
    misc_addVectorsInPlace = &misc_addVectorsInPlace_sse2;
    misc_subtractVectorsInPlace = &misc_subtractVectorsInPlace_sse2;
    misc_addVectorsInPlaceWithMultiplier = &misc_addVectorsInPlaceWithMultiplier_sse2;
    
    misc_addScalarToVectorInPlace = &misc_addScalarToVectorInPlace_sse2;
    misc_setVectorToConstant = & misc_setVectorToConstant_sse2;
    misc_transposeMatrix = &misc_transposeMatrix_sse2;
  } else
#endif
#ifdef COMPILER_SUPPORTS_NEON
  if (i >= MISC_INST_NEON) {
    misc_simd_alignment = 64;
    misc_addAlignedVectorsInPlace = &misc_addAlignedVectorsInPlace_neon;
    misc_subtractAlignedVectorsInPlace = &misc_subtractAlignedVectorsInPlace_neon;
    
    misc_addVectorsInPlace = &misc_addVectorsInPlace_neon;
    misc_subtractVectorsInPlace = &misc_subtractVectorsInPlace_neon;
    misc_addVectorsInPlaceWithMultiplier = &misc_addVectorsInPlaceWithMultiplier_neon;
    
    misc_addScalarToVectorInPlace = & misc_addScalarToVectorInPlace_neon;
    misc_setVectorToConstant = &misc_setVectorToConstant_neon;
    misc_transposeMatrix = &misc_transposeMatrix_neon;
  } else
#endif
  {
    misc_simd_alignment = 0;
    misc_addAlignedVectorsInPlace = &misc_addVectorsInPlace_c;
    misc_subtractAlignedVectorsInPlace = &misc_subtractVectorsInPlace_c;
    
    misc_addVectorsInPlace = &misc_addVectorsInPlace_c;
    misc_subtractVectorsInPlace = &misc_subtractVectorsInPlace_c;
    misc_addVectorsInPlaceWithMultiplier = &misc_addVectorsInPlaceWithMultiplier_c;
    
    misc_addScalarToVectorInPlace = & misc_addScalarToVectorInPlace_c;
    misc_setVectorToConstant = &misc_setVectorToConstant_c;
    misc_transposeMatrix = &misc_transposeMatrix_c;
  }
  
  misc_stat_setSIMDInstructionSet(i);
}

