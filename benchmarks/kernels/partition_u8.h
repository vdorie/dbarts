// u8-width partition kernels for the hot-layer-u8 phase-1 measurement
// (docs/plans/hot-layer-u8.md). Scalar reference plus one SIMD arm (NEON);
// not linked into misc.a, kept local to this benchmark.

#ifndef PARTITION_U8_H
#define PARTITION_U8_H

#include <stddef.h>
#include <stdint.h>

size_t partitionRange_u8_c(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length);
size_t partitionIndices_u8_c(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length);

#if defined(__arm__) || defined(__aarch64__)
size_t partitionRange_u8_neon(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length);
size_t partitionIndices_u8_neon(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length);
#endif

#endif // PARTITION_U8_H
