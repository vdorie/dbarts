#ifndef EXTERNAL_BLOB_IO_H
#define EXTERNAL_BLOB_IO_H

#include "stddef.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
  
// This is pretty rudimentary. It assumes that doubles are IEEE 754 Annex F IEC 60559.
// A byte has to be 8 bits. On the other hand, size_ts are up-cast to 64 bit unsigned integers, 
// and big-endian ints are converted to little-endian first. It should work in most places
// that matter for R. 

struct ext_blob_io;

// open flags are in fcntl.h; permissions in sys/stat.h
ext_blob_io* ext_bio_create(const char* fileName, int openFlag, int permissionsFlag);
void ext_bio_destroy(ext_blob_io* bio);
int ext_bio_initialize(ext_blob_io* bio, const char* fileName, int openFlag, int permissionsFlag);
void ext_bio_invalidate(ext_blob_io* bio);

// write specifies the length before the array is written; writeN just writes N items
// returns 0 on success, otherwise an error code
int ext_bio_writeChars(ext_blob_io* bio, const char* c, size_t length);
int ext_bio_writeNChars(ext_blob_io* bio, const char* c, size_t length);
int ext_bio_writeSizeTypes(ext_blob_io* bio, const size_t* v, size_t length);
int ext_bio_writeNSizeTypes(ext_blob_io* bio, const size_t* v, size_t length);
int ext_bio_writeUnsigned32BitIntegers(ext_blob_io* bio, const uint32_t* v, size_t length);
int ext_bio_writeNUnsigned32BitIntegers(ext_blob_io* bio, const uint32_t* v, size_t length);
int ext_bio_writeDoubles(ext_blob_io* bio, const double* d, size_t length);
int ext_bio_writeNDoubles(ext_blob_io* bio, const double* d, size_t length);

// these return NULL if there is an error and set errno appropriately
char* ext_bio_readChars(ext_blob_io* bio, size_t* length);
char* ext_bio_readNChars(ext_blob_io* bio, size_t length);
size_t* ext_bio_readSizeTypes(ext_blob_io* bio, size_t* length);
size_t* ext_bio_readNSizeTypes(ext_blob_io* bio, size_t length);
uint32_t* ext_bio_readUnsigned32BitIntegers(ext_blob_io* bio, size_t* length);
uint32_t* ext_bio_readNUnsigned32BitIntegers(ext_blob_io* bio, size_t length);
double* ext_bio_readDoubles(ext_blob_io* bio, size_t* length);
double* ext_bio_readNDoubles(ext_blob_io* bio, size_t length);

#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_BLOB_IO_H