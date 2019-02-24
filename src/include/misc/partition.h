#ifndef MISC_PARTITION_H
#define MISC_PARTITION_H

#include <misc/stddef.h>

#include <misc/types.h>

#ifdef __cplusplus
extern "C" {
#endif

extern misc_size_t (*misc_partitionRange)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);
extern misc_size_t (*misc_partitionIndices)(const misc_xint_t* restrict x, misc_xint_t cut, misc_size_t* restrict indices, misc_size_t length);

#ifdef __cplusplus
}
#endif

#endif // MISC_PARTITION_H

