#ifndef MISC_SIMD_H
#define MISC_SIMD_H

#ifdef __cplusplus
extern "C" {
#endif

// call when loading library or else a) segfault and b) no CPU optimized functions will be installed
void misc_simd_init(void);

typedef enum {
  MISC_INST_C = 0,
  MISC_INST_SSE,
  MISC_INST_SSE2,
  MISC_INST_SSE3,
  MISC_INST_SSE4_1,
  MISC_INST_SSE4_2,
  MISC_INST_AVX,
  MISC_INST_AVX2,
  MISC_INST_AVX512F
} misc_simd_instructionLevel;

// THIS IS NOT THREAD SAFE
void misc_simd_setSIMDInstructionSet(misc_simd_instructionLevel i);

misc_simd_instructionLevel misc_simd_getMaxSIMDInstructionSet(void);

#ifdef __cplusplus
}
#endif

#endif // MISC_SIMD_H

