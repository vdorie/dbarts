#ifndef MISC_SIMD_H
#define MISC_SIMD_H

#ifdef __cplusplus
extern "C" {
#endif

// call when loading library or else a) segfault and b) no CPU optimized functions will be installed
void misc_simd_init(void);

typedef enum {
  MISC_INST_C = 0,
  MISC_INST_SSE,      // 1
  MISC_INST_SSE2,     // 2
  MISC_INST_SSE3,     // 3
  MISC_INST_SSSE3,    // 4
  MISC_INST_SSE4_1,   // 5
  MISC_INST_SSE4_2,   // 6
  MISC_INST_AVX,      // 7
  MISC_INST_AVX2,     // 8
  MISC_INST_AVX512F,  // 9
  MISC_INST_AVX512VL, // 10
  MISC_INST_AVX512BW, // 11
  MISC_INST_INVALID
} misc_simd_instructionSet;

// THIS IS NOT THREAD SAFE
void misc_simd_setSIMDInstructionSet(misc_simd_instructionSet i);

misc_simd_instructionSet misc_simd_getMaxSIMDInstructionSet(void);

#ifdef __cplusplus
}
#endif

#endif // MISC_SIMD_H

