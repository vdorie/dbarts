#include "config.h"
#include <misc/partition.h>

#include <stdbool.h>

#include <misc/stddef.h>

#include <misc/intrinsic.h>

#define __USE_NEON__ 1
// from https://stackoverflow.com/questions/11870910/sse-mm-movemask-epi8-equivalent-method-for-arm-neon

/* Input is 16 8 bit unsigned ints, and we want to "movemask", that is 
   "Create mask from the most significant bit of each 8-bit element". 
 */
int vmovemask_u8(uint8x16_t input)
 {
  const int8_t __attribute__ ((aligned (16))) scShift[] = {-7,-6,-5,-4,-3,-2,-1,0,-7,-6,-5,-4,-3,-2,-1,0};
  int8x16_t vshift = vld1q_s8(scShift);
  uint8x16_t vmask = vandq_u8(input, vdupq_n_u8(0x80));
  int32_t out;
   
  vmask = vshlq_u8(vmask, vshift);
  out = vaddv_u8(vget_low_u8(vmask));
  out += (vaddv_u8(vget_high_u8(vmask)) << 8);

  return (int) out;
}

/* Technically, we only need to move masks for 16 bit widths, but the Intel
   intrinsics always operate on 8. Since we are implementing our own we can
   do what we want.
 */

int vmovemask_u16(uint16x8_t input)
{
  const int16_t __attribute__ ((aligned (16))) scShift[] = {-15,-14,-13,-12,-15,-14,-13,-12};
  int16x8_t vshift = vld1q_s16(scShift);
  uint16x8_t vmask = vandq_u16(input, vdupq_n_u16(0x8000));
  // vmask is a packed vector of just 1s if inputs have high bits
  int32_t out;
  
  vmask = vshlq_u16(vmask, vshift);
  // mask now contains 1, 2, 4, 8, 1, 2, 4, 8 if all input are true
  out = vaddv_u16(vget_low_u16(vmask));
  out += (vaddv_u16(vget_high_u16(vmask)) << 4);
  
  return (int) out;
}

#define PARTITION_RANGE 1
size_t misc_partitionRange_neon(const misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#include "partition_body.c"
}

#undef PARTITION_RANGE
#define PARTITION_RANGE 0

size_t misc_partitionIndices_neon(const misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#undef PARTITION_RANGE
#undef __USE_NEON__

