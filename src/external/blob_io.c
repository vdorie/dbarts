#include <ext/blob_io.h>


#include <sys/stat.h> // file permissions
#include <fcntl.h>    // open
#include <unistd.h>   // close, write, sysconf
#include <errno.h>
#include <limits.h>

#ifndef WORDS_BIGENDIAN
#define XOR_SWAP(_X_, _Y_) { (_X_) ^= (_Y_); (_Y_) ^= (_X_); (_X_) ^= (_Y_); }
static inline void swapEndiannessFor32BitInt(uint32_t* u)
{
  char* c = (char*) u;
  XOR_SWAP(c[0], c[3]);
  XOR_SWAP(c[1], c[2]);
}
static inline void swapEndiannessFor64BitInt(uint64_t* u)
{
  char* c = (char*) u;
  XOR_SWAP(c[0], c[7]);
  XOR_SWAP(c[1], c[6]);
  XOR_SWAP(c[2], c[5]);
  XOR_SWAP(c[3], c[4]);
}

static void swapEndiannessFor32BitInts(uint64_t* u, size_t length);
static void swapEndiannessFor64BitInts(uint64_t* u, size_t length);
#endif

typedef struct
{
  int fileDescriptor;
  char* buffer;
  size_t bufferLength;
  
} ext_blob_io;

static int writeSizeType(ext_blob_io* bio, size_t s);
static size_t fillBufferFromSizeTypes(ext_blob_io* restrict bio, const size_t* restrict v, size_t length);

ext_blob_io* ext_bio_create(const char* fileName, int openFlag, int permissionsFlag)
{
  ext_blob_io* result = (ext_blob_io*) malloc(sizeof(ext_blob_io)); // sets ernno = ENOMEM
  if (result == NULL) return NULL;
  
  int errorCode = ext_bio_initialize(result, fileName, openFlag, permissionsFlag);
  if (errorCode != 0) {
    free(result);
    errno = errorCode;
    return NULL;
  }
  
  return result;
}

void ext_bio_destroy(ext_blob_io* bio)
{
  if (bio == NULL) return;
  
  ext_bio_invalidate(bio);
  free(bio);
}

int ext_bio_initialize(ext_blob_io* bio, const char* fileName, int openFlag, int permissionsFlag)
{
  if (bio == NULL) return EFAULT;
  
  bio->buffer = NULL;
  bio->bufferLength = 0;
  
  bio->fileDescriptor = open(fileName, openFlag, permissionsFlag);
  if (bio->fileDescriptor == -1) return errno;
  
  errno = 0;
  long pageSize = sysconf(_SC_PAGE_SIZE);
  if (pageSize <= 0 || errno != 0) pageSize = PAGE_SIZE; // is set when compiled, but better than nothing I guess
  bio->bufferLength = pageSize;
    
  bio->buffer = (char*) malloc(bufferLength);
  if (bio->buffer == NULL) {
    close(bio->fileDescriptor);
    bio->fileDescriptor = -1;
    bio->bufferLength = 0;
    return ENOMEM;
  }
  
  return 0;
}

void ext_bio_invalidate(ext_blob_io* bio)
{
  if (bio == NULL) return;
  
  if (bio->fileDescriptor != -1) {
    close(bio->fileDescriptor);
    bio->fileDescriptor = -1;
  }
  
  if (bio->buffer != NULL) {
    free(bio->buffer);
    bio->buffer = NULL;
  }
  
  bio->bufferLength = 0;
}

static int writeSizeType(ext_blob_io* bio, size_t s)
{
  uint64_t u = (uint64_t) s;
#ifndef WORDS_BIGENDIAN
  swapEndiannessForInt(&u);
#endif
  
  ssize_t bytesWritten = write(bio->fileDescriptor, &u, sizeof(uint64_t));
  if (bytesWritten == 0) return EIO;
  if (bytesWritten < 0) return errno;
  
  return 0;
}

int ext_bio_writeChars(ext_blob_io* bio, const char* c, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  int errorCode = writeSizeType(bio, length);
  if (errorCode != 0) return result;
  
  return ext_bio_writeNChars(bio, c, length);
}

int ext_bio_writeNChars(ext_blob_io* bio, const char* c, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  size_t totalBytesWritten = 0;
  while (totalBytesWritten < length) {
    ssize_t bytesWritten = write(bio->fileDescriptor, c + totalBytesWritten, length - totalBytesWritten);
    if (bytesWritten == 0) return EIO;
    if (bytesWritten < 0) return errno;
    totalBytesWritten += bytesWritten;
  }
  
  return 0;
}

int ext_bio_writeSizeTypes(ext_blob_io* bio, const size_t* v, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  int errorCode = writeSizeType(bio, length);
  if (errorCode != 0) return result;
  
  return ext_bio_writeNSizeTypes(bio, v, length);
}

int ext_bio_writeNSizeTypes(ext_blob_io* bio, const size_t* v, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  size_t totalItemsWritten = 0;
  while (totalItemsWritten < length) {
    size_t itemsWritten = fillBufferFromSizeTypes(bio, v + totalItemsWritten, length - totalItemsWritten);
    
    size_t totalBytesWritten = 0;
    size_t bytesToWrite = itemsWritten * sizeof(uint64_t);
    
    while (totalBytesWritten < bytesToWrite) {
      ssize_t bytesWritten = write(bio->fileDescriptor, bio->buffer, bytesToWrite - totalBytesWritten);
      if (bytesWritten == 0) return EIO;
      if (bytesWritten < 0) return errno;
      totalBytesWritten += bytesWritten;
    }
  }
  
  return 0;
}

int ext_bio_writeUnsigned32BitIntegers(ext_blob_io* bio, const uint32_t* v, size_t length)
{
  uint64_t lengthInt = (uint64_t) length;
#ifndef WORDS_BIGENDIAN
  swapEndiannessForInt(&lengthInt);
#endif
  
  ssize_t bytesWritten = write(bio->fileDescriptor, &lengthInt, sizeof(uint64_t));
  if (bytesWritten == 0) return EIO;
  if (bytesWritten < 0) return errno;
  
  return ext_bio_writeNUnsigned32BitIntegers(bio, v, length);
}

int ext_bio_writeNUnsigned32BitIntegers(ext_blob_io* bio, const uint32_t* v, size_t length)
{
  size_t totalItemsWritten = 0;
  while (totalItemsWritten < length) {
    size_t itemsWritten = fillBufferFrom32BitIntegers(bio, v + totalItemsWritten, length - totalItemsWritten);
    
    size_t totalBytesWritten = 0;
    size_t bytesToWrite = itemsWritten * sizeof(uint32_t);
    
    while (totalBytesWritten < bytesToWrite) {
      ssize_t bytesWritten = write(bio->fileDescriptor, bio->buffer, bytesToWrite - totalBytesWritten);
      if (bytesWritten == 0) return EIO;
      if (bytesWritten < 0) return errno;
      totalBytesWritten += bytesWritten;
    }
  }
  
  return 0;
}

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

// upcast to 64_bit ints
static size_t fillBufferFromSizeTypes(ext_blob_io* restrict bio, const size_t* restrict v, size_t length)
{
  if (length == 0) return 0;
  
  uint64_t* restrict buffer = (uint64_t* restrict) bio->buffer;
  size_t bufferLength = bio->bufferLength / sizeof(uint64_t);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
  if (sizeof(size_t) == sizeof(uint64_t)) {
    memcpy(buffer, v, fillLength * sizeof(uint64_t));
#ifndef WORDS_BIGENDIAN
   swapEndiannessFor64BitInts(buffer, billLength);
#endif
    return fillLength;
  }

  size_t lengthMod5 = fillLength % 5;
  size_t i = 0;
  for ( ; i < lengthMod5; ++i) buffer[i] = (uint64_t) v[i];
  
  for ( ; i < fillLength; i += 5) {
    buffer[    i]     = (uint64_t) v[i];
    buffer[i + 1] = (uint64_t) v[i + 1];
    buffer[i + 2] = (uint64_t) v[i + 2];
    buffer[i + 3] = (uint64_t) v[i + 3];
    buffer[i + 4] = (uint64_t) v[i + 4];
  }
  
#ifndef WORDS_BIGENDIAN
   swapEndiannessFor64BitInts(buffer, fillLength);
#endif
  
  return fillLength;
}

static size_t fillBufferFromUnsigned32BitIntegers(ext_blob_io* restrict bio, const uint32_t* restrict v, size_t length)
{
  if (length == 0) return 0;
  
  uint32_t* restrict buffer = (uint32_t* restrict) bio->buffer;
  size_t bufferLength = bio->bufferLength / sizeof(uint32_t);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
  memcpy(buffer, v, fillLength * sizeof(uint32_t));
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor32BitInts(buffer, billLength);
#endif
   
  return fillLength;
}

#ifndef WORDS_BIGENDIAN
static void swapEndiannessFor32BitInts(uint32_t* v, size_t length)
{
  size_t lengthMod5 = length % 5;
  size_t i = 0;
  for ( ; i < lengthMod5; ++i) swapEndiannessFor32BitInt(v + i);
  
  for ( ; i < length; i += 5) {
    swapEndiannessFor32BitInt(v + i);
    swapEndiannessFor32BitInt(v + i + 1);
    swapEndiannessFor32BitInt(v + i + 2);
    swapEndiannessFor32BitInt(v + i + 3);
    swapEndiannessFor32BitInt(v + i + 4);
  }
}

static void swapEndiannessFor64BitInts(uint64_t* v, size_t length)
{
  size_t lengthMod5 = length % 5;
  size_t i = 0;
  for ( ; i < lengthMod5; ++i) swapEndiannessFor64BitInt(v + i);
  
  for ( ; i < length; i += 5) {
    swapEndiannessFor64BitInt(v + i);
    swapEndiannessFor64BitInt(v + i + 1);
    swapEndiannessFor64BitInt(v + i + 2);
    swapEndiannessFor64BitInt(v + i + 3);
    swapEndiannessFor64BitInt(v + i + 4);
  }
}
#endif