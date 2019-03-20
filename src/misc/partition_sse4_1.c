#include "config.h"
#include <misc/partition.h>

#include <stdbool.h>

#include <misc/stddef.h>

#include <misc/intrinsic.h>

#define __USE_SSE4_1__ 1

#define PARTITION_RANGE 1
size_t misc_partitionRange_sse4_1(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#include "partition_body.c"
}

#undef PARTITION_RANGE
#define PARTITION_RANGE 0

size_t misc_partitionIndices_sse4_1(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#undef PARTITION_RANGE
#undef __USE_SSE4_1__

