#include "config.h"
#include <misc/partition.h>

#include <stdbool.h>

#include <misc/stddef.h>

#include <misc/intrinsic.h>

#define __USE_NEON__ 1
// from https://stackoverflow.com/questions/11870910/sse-mm-movemask-epi8-equivalent-method-for-arm-neon
int vmovemask_u8(uint16x8_t input)
{
  const uint8_t __attribute__ ((aligned (16))) ucShift[] = {-7,-6,-5,-4,-3,-2,-1,0,-7,-6,-5,-4,-3,-2,-1,0};
  uint8x16_t vshift = vld1q_u8(ucShift);
  uint8x16_t vmask = vandq_u8(input, vdupq_n_u8(0x80));
  uint32_t out;
  
  vmask = vshlq_u8(vmask, vshift);
  out = vaddv_u8(vget_low_u8(vmask));
  out += (vaddv_u8(vget_high_u8(vmask)) << 8);
  
  return out;
}

#define PARTITION_RANGE 1
size_t misc_partitionRange_neon(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#include "partition_body.c"
}

#undef PARTITION_RANGE
#define PARTITION_RANGE 0

size_t misc_partitionIndices_neon(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#undef PARTITION_RANGE
#undef __USE_NEON__

