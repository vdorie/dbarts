#ifndef EXTERNAL_BINARY_IO_H
#define EXTERNAL_BINARY_IO_H

#include "stddef.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
  
// This is pretty rudimentary. It assumes that doubles are IEEE 754 Annex F IEC 60559.
// A byte has to be 8 bits. On the other hand, size_ts are up-cast to 64 bit unsigned integers, 
// and big-endian ints are converted to little-endian first. It should work in most places
// that matter for R. 

typedef struct ext_binaryIO
{
  int fileDescriptor;
  char* buffer;
  size_t bufferLength;
  
} ext_binaryIO;

// open flags are in fcntl.h; permissions in sys/stat.h
ext_binaryIO* ext_bio_create(const char* fileName, int openFlag, int permissionsFlag);
void ext_bio_destroy(ext_binaryIO* bio);
int ext_bio_initialize(ext_binaryIO* bio, const char* fileName, int openFlag, int permissionsFlag);
void ext_bio_invalidate(ext_binaryIO* bio);

// write specifies the length before the array is written; writeN just writes N items
// returns 0 on success, otherwise an error code
int ext_bio_writeChar(ext_binaryIO* bio, char c);
int ext_bio_writeChars(ext_binaryIO* bio, const char* c, size_t length);
int ext_bio_writeNChars(ext_binaryIO* bio, const char* c, size_t length);
int ext_bio_writeSizeType(ext_binaryIO* bio, size_t s);
int ext_bio_writeSizeTypes(ext_binaryIO* bio, const size_t* s, size_t length);
int ext_bio_writeNSizeTypes(ext_binaryIO* bio, const size_t* s, size_t length);
int ext_bio_writeUnsigned32BitInteger(ext_binaryIO* bio, uint32_t u);
int ext_bio_writeUnsigned32BitIntegers(ext_binaryIO* bio, const uint32_t* u, size_t length);
int ext_bio_writeNUnsigned32BitIntegers(ext_binaryIO* bio, const uint32_t* u, size_t length);
int ext_bio_writeUnsigned64BitInteger(ext_binaryIO* bio, uint64_t u);
int ext_bio_writeDouble(ext_binaryIO* bio, double d);
int ext_bio_writeDoubles(ext_binaryIO* bio, const double* d, size_t length);
int ext_bio_writeNDoubles(ext_binaryIO* bio, const double* d, size_t length);

// for anything that returns a pointer, if it fails it returns NULL and sets errno
int ext_bio_readChar(ext_binaryIO* bio, char* c);
char* ext_bio_readChars(ext_binaryIO* bio, size_t* length);
int ext_bio_readNChars(ext_binaryIO* bio, char* c, size_t length);
int ext_bio_readSizeType(ext_binaryIO* bio, size_t* s);
size_t* ext_bio_readSizeTypes(ext_binaryIO* bio, size_t* length);
int ext_bio_readNSizeTypes(ext_binaryIO* bio, size_t* s, size_t length);
int ext_bio_readUnsigned32BitInteger(ext_binaryIO* bio, uint32_t* u);
uint32_t* ext_bio_readUnsigned32BitIntegers(ext_binaryIO* bio, size_t* length);
int ext_bio_readNUnsigned32BitIntegers(ext_binaryIO* bio, uint32_t* u, size_t length);
int ext_bio_readUnsigned64BitInteger(ext_binaryIO* bio, uint64_t *u);
int ext_bio_readDouble(ext_binaryIO* bio, double* d);
double* ext_bio_readDoubles(ext_binaryIO* bio, size_t* length);
int ext_bio_readNDoubles(ext_binaryIO* bio, double* d, size_t length);

#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_BINARY_IO_H
