#include "config.h"
#include <misc/partition.h>

#include <stdbool.h>

#include <misc/stddef.h>

#include <misc/intrinsic.h>


#include <external/io.h>

misc_size_t (*misc_partitionRange)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length) = 0;
misc_size_t (*misc_partitionIndices)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length) = 0;

#ifdef __AVX2__
#  define __USE_AVX2__ 1

#  define PARTITION_RANGE 1
size_t misc_partitionRange_avx2(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#  undef PARTITION_RANGE
#  define PARTITION_RANGE 0

size_t misc_partitionIndices_avx2(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#  undef PARTITION_RANGE
#  undef __USE_AVX2__
#endif

#ifdef __SSE4_1__
#  define __USE_SSE4_1__ 1

#  define PARTITION_RANGE 1
size_t misc_partitionRange_sse4_1(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#  undef PARTITION_RANGE
#  define PARTITION_RANGE 0

size_t misc_partitionIndices_sse4_1(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#  undef PARTITION_RANGE
#  undef __USE_SSE4_1__
#endif

#ifdef __SSE2__
#  define __USE_SSE2__ 1

#  define PARTITION_RANGE 1
size_t misc_partitionRange_sse2(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#  undef PARTITION_RANGE
#  define PARTITION_RANGE 0

size_t misc_partitionIndices_sse2(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#  include "partition_body.c"
}

#  undef PARTITION_RANGE
#  undef __USE_SSE2__
#endif

#define PARTITION_RANGE 1
size_t misc_partitionRange_c(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#include "partition_body.c"
}

#undef PARTITION_RANGE
#define PARTITION_RANGE 0

size_t misc_partitionIndices_c(misc_xint_t* restrict x, misc_xint_t cut, size_t* restrict indices, size_t length)
{
#include "partition_body.c"
}

#undef PARTITION_RANGE

