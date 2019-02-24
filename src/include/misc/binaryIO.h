#ifndef MISC_BINARY_IO_H
#define MISC_BINARY_IO_H

#include <misc/stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
  
// This is pretty rudimentary. It assumes that doubles are IEEE 754 Annex F IEC 60559.
// A byte has to be 8 bits. On the other hand, misc_size_ts are up-cast to 64 bit unsigned integers, 
// and big-endian ints are converted to little-endian first. It should work in most places
// that matter for R. 

typedef struct misc_binaryIO
{
  int fileDescriptor;
  void* buffer;
  misc_size_t bufferLength;
  
} misc_binaryIO;

// open flags are in fcntl.h; permissions in sys/stat.h
misc_binaryIO* misc_bio_create(const char* fileName, int openFlag, int permissionsFlag);
void misc_bio_destroy(misc_binaryIO* bio);
int misc_bio_initialize(misc_binaryIO* bio, const char* fileName, int openFlag, int permissionsFlag);
void misc_bio_invalidate(misc_binaryIO* bio);

// write specifies the length before the array is written; writeN just writes N items
// returns 0 on success, otherwise an error code
int misc_bio_writeChar(misc_binaryIO* bio, char c);
int misc_bio_writeChars(misc_binaryIO* bio, const char* c, misc_size_t length);
int misc_bio_writeNChars(misc_binaryIO* bio, const char* c, misc_size_t length);
int misc_bio_writeSizeType(misc_binaryIO* bio, misc_size_t s);
int misc_bio_writeSizeTypes(misc_binaryIO* bio, const misc_size_t* s, misc_size_t length);
int misc_bio_writeNSizeTypes(misc_binaryIO* bio, const misc_size_t* s, misc_size_t length);
int misc_bio_writeUnsigned32BitInteger(misc_binaryIO* bio, uint32_t u);
int misc_bio_writeUnsigned32BitIntegers(misc_binaryIO* bio, const uint32_t* u, misc_size_t length);
int misc_bio_writeNUnsigned32BitIntegers(misc_binaryIO* bio, const uint32_t* u, misc_size_t length);
int misc_bio_writeUnsigned64BitInteger(misc_binaryIO* bio, uint64_t u);
int misc_bio_writeDouble(misc_binaryIO* bio, double d);
int misc_bio_writeDoubles(misc_binaryIO* bio, const double* d, misc_size_t length);
int misc_bio_writeNDoubles(misc_binaryIO* bio, const double* d, misc_size_t length);

// takes ints but writes them as 64 bit unsigned
int misc_bio_writeNInts(misc_binaryIO* bio, const int* i, misc_size_t length);

// for anything that returns a pointer, if it fails it returns NULL and sets errno
int misc_bio_readChar(misc_binaryIO* bio, char* c);
char* misc_bio_readChars(misc_binaryIO* bio, misc_size_t* length);
int misc_bio_readNChars(misc_binaryIO* bio, char* c, misc_size_t length);
int misc_bio_readSizeType(misc_binaryIO* bio, misc_size_t* s);
size_t* misc_bio_readSizeTypes(misc_binaryIO* bio, misc_size_t* length);
int misc_bio_readNSizeTypes(misc_binaryIO* bio, misc_size_t* s, misc_size_t length);
int misc_bio_readUnsigned32BitInteger(misc_binaryIO* bio, uint32_t* u);
uint32_t* misc_bio_readUnsigned32BitIntegers(misc_binaryIO* bio, misc_size_t* length);
int misc_bio_readNUnsigned32BitIntegers(misc_binaryIO* bio, uint32_t* u, misc_size_t length);
int misc_bio_readUnsigned64BitInteger(misc_binaryIO* bio, uint64_t *u);
int misc_bio_readDouble(misc_binaryIO* bio, double* d);
double* misc_bio_readDoubles(misc_binaryIO* bio, misc_size_t* length);
int misc_bio_readNDoubles(misc_binaryIO* bio, double* d, misc_size_t length);

// reads in ints from 64 bit unsigned - can result in loss of precision
int misc_bio_readNInts(misc_binaryIO* bio, int* i, misc_size_t length);

#ifdef __cplusplus
}
#endif

#endif // MISC_BINARY_IO_H

