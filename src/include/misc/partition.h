#ifndef MISC_PARTITION_H
#define MISC_PARTITION_H

#include <stdint.h>

#include <misc/stddef.h>

#include <misc/types.h>

#ifdef __cplusplus
extern "C" {
#endif

extern misc_size_t (*misc_partitionRange)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern misc_size_t (*misc_partitionIndices)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);

// The rank-bitmap sparse-column partition (docs/design/sparse-columns.md):
// code(i) is zeroCode when bit i of bits is clear, otherwise
// nzCodes[wordRanks[i / 64] + popcount of the lower bits of word i / 64].
// Same contract as misc_partitionIndices; scalar only, so no dispatch
// pointer yet.
misc_size_t misc_partitionIndicesSparse(const uint64_t* restrict bits, const uint32_t* restrict wordRanks, const misc_xint_t* restrict nzCodes, misc_xint_t zeroCode, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);

#ifdef __cplusplus
}
#endif

#endif // MISC_PARTITION_H

