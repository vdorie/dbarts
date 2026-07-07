#include "partition_u8.h"

#include <stdbool.h>

// One Hoare two-pointer step; gather selects x[idx[i]] (misc_partitionIndices,
// idx already a permutation) over x[i] (misc_partitionRange, idx seeded by
// the caller). Returns false once lh/rh meet, without swapping.
static inline bool scalarStep_u8(const uint8_t* restrict x, uint8_t cut, size_t* restrict idx, size_t* lh, size_t* rh, bool gather)
{
  while ((gather ? x[idx[*lh]] : x[*lh]) <= cut && *lh < *rh) ++*lh;
  while ((gather ? x[idx[*rh]] : x[*rh])  > cut && *lh < *rh) --*rh;
  if (*lh >= *rh) return false;

  if (gather) { size_t t = idx[*rh]; idx[*rh] = idx[*lh]; idx[*lh] = t; }
  else        { idx[*rh] = *lh; idx[*lh] = *rh; }
  ++*lh; --*rh;
  return true;
}

static inline size_t finish_u8(const uint8_t* restrict x, uint8_t cut, const size_t* restrict idx, size_t lh, bool gather)
{
  return (gather ? x[idx[lh]] : x[lh]) <= cut ? lh + 1 : lh;
}

static size_t partitionCore_u8_c(const uint8_t* restrict x, uint8_t cut, size_t* restrict idx, size_t length, bool gather)
{
  if (length == 0) return 0;
  size_t lh = 0, rh = length - 1;
  while (scalarStep_u8(x, cut, idx, &lh, &rh, gather)) { }
  return finish_u8(x, cut, idx, lh, gather);
}

size_t partitionRange_u8_c(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length)
{
  for (size_t i = 0; i < length; ++i) indices[i] = i;
  return partitionCore_u8_c(x, cut, indices, length, false);
}

size_t partitionIndices_u8_c(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length)
{
  return partitionCore_u8_c(x, cut, indices, length, true);
}

#if defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>

#define U8_BLOCK 16

static inline uint8x16_t gather16_u8(const uint8_t* restrict x, const size_t* restrict idx, size_t pos)
{
  uint8_t buf[U8_BLOCK];
  for (size_t k = 0; k < U8_BLOCK; ++k) buf[k] = x[idx[pos + k]];
  return vld1q_u8(buf);
}

// Same two-pointer step as partitionCore_u8_c, but a block of 16 is skipped
// with one vector reduce (vmaxvq/vminvq) whenever it is already entirely on
// the correct side; a mismatch drops one scalar step before re-probing.
static size_t partitionCore_u8_neon(const uint8_t* restrict x, uint8_t cut, size_t* restrict idx, size_t length, bool gather)
{
  if (length == 0) return 0;
  size_t lh = 0, rh = length - 1;

  while (lh + 2 * U8_BLOCK <= rh) {
    uint8x16_t lv = gather ? gather16_u8(x, idx, lh) : vld1q_u8(x + lh);
    if (vmaxvq_u8(lv) <= cut) { lh += U8_BLOCK; continue; }

    uint8x16_t rv = gather ? gather16_u8(x, idx, rh - (U8_BLOCK - 1)) : vld1q_u8(x + rh - (U8_BLOCK - 1));
    if (vminvq_u8(rv) > cut) { rh -= U8_BLOCK; continue; }

    if (!scalarStep_u8(x, cut, idx, &lh, &rh, gather)) return finish_u8(x, cut, idx, lh, gather);
  }

  while (scalarStep_u8(x, cut, idx, &lh, &rh, gather)) { }
  return finish_u8(x, cut, idx, lh, gather);
}

size_t partitionRange_u8_neon(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length)
{
  for (size_t i = 0; i < length; ++i) indices[i] = i;
  return partitionCore_u8_neon(x, cut, indices, length, false);
}

size_t partitionIndices_u8_neon(const uint8_t* restrict x, uint8_t cut, size_t* restrict indices, size_t length)
{
  return partitionCore_u8_neon(x, cut, indices, length, true);
}

#undef U8_BLOCK
#endif // __arm__ || __aarch64__
