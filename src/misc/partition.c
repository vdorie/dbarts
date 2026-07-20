#include "config.h"
#include <misc/partition.h>

#include <stdbool.h>

#include <misc/stddef.h>

misc_size_t (*misc_partitionRange)(const misc_xint_t* restrict x, misc_xint_t cut, misc_index_t* restrict indices, misc_size_t length) = 0;
misc_size_t (*misc_partitionIndices)(const misc_xint_t* restrict x, misc_xint_t cut, misc_index_t* restrict indices, misc_size_t length) = 0;

#define PARTITION_RANGE 1
size_t misc_partitionRange_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_index_t* restrict indices, size_t length)
{
#include "partition_body.c"
}

#undef PARTITION_RANGE
#define PARTITION_RANGE 0

size_t misc_partitionIndices_c(const misc_xint_t* restrict x, misc_xint_t cut, misc_index_t* restrict indices, size_t length)
{
#include "partition_body.c"
}

#undef PARTITION_RANGE

#if defined(__GNUC__) || defined(__clang__)
#  define popcount64(_X_) ((unsigned int) __builtin_popcountll(_X_))
#else
static inline unsigned int popcount64(uint64_t v) {
  v = v - ((v >> 1) & UINT64_C(0x5555555555555555));
  v = (v & UINT64_C(0x3333333333333333)) + ((v >> 2) & UINT64_C(0x3333333333333333));
  v = (v + (v >> 4)) & UINT64_C(0x0F0F0F0F0F0F0F0F);
  return (unsigned int) ((v * UINT64_C(0x0101010101010101)) >> 56);
}
#endif

static inline misc_xint_t sparseCodeAt(const uint64_t* restrict bits, const uint32_t* restrict wordRanks, const misc_xint_t* restrict nzCodes, misc_xint_t zeroCode, size_t i)
{
  uint64_t word = bits[i >> 6];
  uint64_t bit = UINT64_C(1) << (i & 63u);
  if ((word & bit) == 0) return zeroCode;
  return nzCodes[wordRanks[i >> 6] + popcount64(word & (bit - 1u))];
}

// the scalar two-pointer loop of partition_body.c over rank-bitmap access
size_t misc_partitionIndicesSparse(const uint64_t* restrict bits, const uint32_t* restrict wordRanks, const misc_xint_t* restrict nzCodes, misc_xint_t zeroCode, misc_xint_t cut, misc_index_t* restrict indices, size_t length)
{
  if (length == 0) return 0;

  size_t lh = 0, rh = length - 1;

  while (true) {
    while (sparseCodeAt(bits, wordRanks, nzCodes, zeroCode, indices[lh]) <= cut && lh < rh) ++lh;
    while (sparseCodeAt(bits, wordRanks, nzCodes, zeroCode, indices[rh])  > cut && lh < rh) --rh;

    if (lh >= rh) break;

    misc_index_t temp = indices[rh];
    indices[rh] = indices[lh];
    indices[lh] = temp;

    ++lh;
    --rh;
  }

  return sparseCodeAt(bits, wordRanks, nzCodes, zeroCode, indices[lh]) <= cut ? lh + 1 : lh;
}

