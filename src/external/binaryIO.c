#include <external/binaryIO.h>

#include <sys/stat.h> // file permissions
#include <fcntl.h>    // open
#include <unistd.h>   // close, write, sysconf
#include <stdlib.h>   // malloc, posix_memalign
#ifdef HAVE_MALLOC_H
#include <malloc.h>   // __mingw_aligned_malloc
#endif
#include <string.h>   // memcpy
#include <errno.h>
#include <limits.h>

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

#ifndef EOVERFLOW // needed for mingw
#define EOVERFLOW EFBIG
#endif

#ifndef WORDS_BIGENDIAN
#define XOR_SWAP(_X_, _Y_) { (_X_) ^= (_Y_); (_Y_) ^= (_X_); (_X_) ^= (_Y_); }

static inline void swapEndiannessFor4ByteWord(char* c)
{
  XOR_SWAP(c[0], c[3]);
  XOR_SWAP(c[1], c[2]);
}

static inline void swapEndiannessFor8ByteWord(char* c)
{
  XOR_SWAP(c[0], c[7]);
  XOR_SWAP(c[1], c[6]);
  XOR_SWAP(c[2], c[5]);
  XOR_SWAP(c[3], c[4]);
}

static void swapEndiannessFor4ByteWords(char* c, size_t length);
static void swapEndiannessFor8ByteWords(char* c, size_t length);
#endif // WORDS_BIGENDIAN

static size_t fillBufferFromSizeTypes(ext_binaryIO* restrict bio, const size_t* restrict v, size_t length);
static size_t fillSizeTypesFromBuffer(ext_binaryIO* restrict bio, size_t* restrict v, size_t length);
static size_t fillBufferFromDoubles(ext_binaryIO* restrict bio, const double* restrict v, size_t length);
static size_t fillDoublesFromBuffer(ext_binaryIO* restrict bio, double* restrict v, size_t length);  
static size_t fillBufferFromUnsigned32BitIntegers(ext_binaryIO* restrict bio, const uint32_t* restrict v, size_t length);
static size_t fillUnsigned32BitIntegersFromBuffer(ext_binaryIO* restrict bio, uint32_t* restrict v, size_t length);

ext_binaryIO* ext_bio_create(const char* fileName, int openFlag, int permissionsFlag)
{
  ext_binaryIO* result = (ext_binaryIO*) malloc(sizeof(ext_binaryIO)); // sets ernno = ENOMEM
  if (result == NULL) return NULL;
  
  int errorCode = ext_bio_initialize(result, fileName, openFlag, permissionsFlag);
  if (errorCode != 0) {
    free(result);
    errno = errorCode;
    return NULL;
  }
  
  return result;
}

void ext_bio_destroy(ext_binaryIO* bio)
{
  if (bio == NULL) return;
  
  ext_bio_invalidate(bio);
  free(bio);
}

int ext_bio_initialize(ext_binaryIO* bio, const char* fileName, int openFlag, int permissionsFlag)
{
  if (bio == NULL) return EFAULT;
  
  bio->buffer = NULL;
  bio->bufferLength = 0;
  
  bio->fileDescriptor = open(fileName, openFlag, permissionsFlag);
  if (bio->fileDescriptor == -1) return errno;
  
  errno = 0;
#ifdef _WIN32
  SYSTEM_INFO si;
  GetSystemInfo(&si);

  long pageSize = si.dwPageSize;
#elif defined(_SC_PAGE_SIZE)
  long pageSize = sysconf(_SC_PAGE_SIZE);
#elif defined(_SC_PAGESIZE)
  long pageSize = sysconf(_SC_PAGESIZE)
#else
  long pageSize = 4096;
#endif

  if (pageSize <= 0 || errno != 0) pageSize = 4096; // sure, why not?
  bio->bufferLength = (size_t) pageSize;
  
  // practically needs to hold at least one 64 bit int
  while (bio->bufferLength < sizeof(uint64_t)) bio->bufferLength <<= 1;
  
  size_t alignment = sizeof(uint64_t);
  if (alignment % sizeof(void*) != 0) alignment *= sizeof(void*);
  
#ifdef __MINGW32__
  bio->buffer = __mingw_aligned_malloc(bio->bufferLength, alignment);
  if (bio->buffer == NULL) {
    int errorCode = ENOMEM;
#else
  int errorCode = posix_memalign(&bio->buffer, alignment, bio->bufferLength);
  if (errorCode != 0) {
    if (bio->buffer != NULL) free(bio->buffer);
#endif
    close(bio->fileDescriptor);
    bio->fileDescriptor = -1;
    bio->bufferLength = 0;
    return errorCode;
  }
  
  return 0;
}

void ext_bio_invalidate(ext_binaryIO* bio)
{
  if (bio == NULL) return;
  
  if (bio->fileDescriptor != -1) {
    close(bio->fileDescriptor);
    bio->fileDescriptor = -1;
  }
  
  if (bio->buffer != NULL) {
#ifdef __MINGW32__
    __mingw_aligned_free(bio->buffer);
#else
    free(bio->buffer);
#endif
    bio->buffer = NULL;
  }
  
  bio->bufferLength = 0;
}

int ext_bio_writeChar(ext_binaryIO* bio, char c)
{
  if (bio == NULL) return EFAULT;
  
  ssize_t bytesWritten = write(bio->fileDescriptor, &c, sizeof(char));
  
  if (bytesWritten == 0) return EIO;
  if (bytesWritten < 0) return errno;
  
  return 0;
}

int ext_bio_writeChars(ext_binaryIO* bio, const char* c, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  int errorCode = ext_bio_writeSizeType(bio, length);
  if (errorCode != 0) return errorCode;
  
  return ext_bio_writeNChars(bio, c, length);
}

int ext_bio_writeNChars(ext_binaryIO* bio, const char* c, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  size_t totalBytesWritten = 0;
  while (totalBytesWritten < length) {
    ssize_t bytesWritten = write(bio->fileDescriptor, c + totalBytesWritten, length - totalBytesWritten);
    if (bytesWritten == 0) return EIO;
    if (bytesWritten < 0) return errno;
    totalBytesWritten += (size_t) bytesWritten;
  }
  
  return 0;
}

int ext_bio_writeSizeType(ext_binaryIO* bio, size_t s)
{
  uint64_t u = (uint64_t) s;
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor8ByteWord((char*) &u);
#endif
  
  ssize_t bytesWritten = write(bio->fileDescriptor, &u, sizeof(uint64_t));
  if (bytesWritten == 0) return EIO;
  if (bytesWritten < 0) return errno;
  
  return 0;
}

int ext_bio_writeSizeTypes(ext_binaryIO* bio, const size_t* v, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  int errorCode = ext_bio_writeSizeType(bio, length);
  if (errorCode != 0) return errorCode;
  
  return ext_bio_writeNSizeTypes(bio, v, length);
}

int ext_bio_writeNSizeTypes(ext_binaryIO* bio, const size_t* v, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  size_t totalItemsWritten = 0;
  while (totalItemsWritten < length) {
    size_t itemsWritten = fillBufferFromSizeTypes(bio, v + totalItemsWritten, length - totalItemsWritten);
    
    size_t totalBytesWritten = 0;
    size_t bytesToWrite = itemsWritten * sizeof(uint64_t);
    
    while (totalBytesWritten < bytesToWrite) {
      ssize_t bytesWritten = write(bio->fileDescriptor, (char*) bio->buffer, bytesToWrite - totalBytesWritten);
      if (bytesWritten == 0) return EIO;
      if (bytesWritten < 0) return errno;
      totalBytesWritten += (size_t) bytesWritten;
    }
    
    totalItemsWritten += itemsWritten;
  }
  
  return 0;
}

int ext_bio_writeUnsigned32BitInteger(ext_binaryIO* bio, uint32_t u)
{
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor4ByteWord((char*) &u);
#endif
  
  ssize_t bytesWritten = write(bio->fileDescriptor, &u, sizeof(uint32_t));
  if (bytesWritten == 0) return EIO;
  if (bytesWritten < 0) return errno;
  
  return 0;
}

int ext_bio_writeUnsigned32BitIntegers(ext_binaryIO* bio, const uint32_t* u, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  int errorCode = ext_bio_writeSizeType(bio, length);
  if (errorCode != 0) return errorCode;
  
  return ext_bio_writeNUnsigned32BitIntegers(bio, u, length);
}

int ext_bio_writeNUnsigned32BitIntegers(ext_binaryIO* bio, const uint32_t* u, size_t length)
{
  size_t totalItemsWritten = 0;
  while (totalItemsWritten < length) {
    size_t itemsWritten = fillBufferFromUnsigned32BitIntegers(bio, u + totalItemsWritten, length - totalItemsWritten);
    
    size_t totalBytesWritten = 0;
    size_t bytesToWrite = itemsWritten * sizeof(uint32_t);
    
    while (totalBytesWritten < bytesToWrite) {
      ssize_t bytesWritten = write(bio->fileDescriptor, (char*) bio->buffer, bytesToWrite - totalBytesWritten);
      if (bytesWritten == 0) return EIO;
      if (bytesWritten < 0) return errno;
      totalBytesWritten += (size_t) bytesWritten;
    }
    
    totalItemsWritten += itemsWritten;
  }
  
  return 0;
}

int ext_bio_writeUnsigned64BitInteger(ext_binaryIO* bio, uint64_t u)
{
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor8ByteWord((char*) &u);
#endif
  
  ssize_t bytesWritten = write(bio->fileDescriptor, &u, sizeof(uint64_t));
  if (bytesWritten == 0) return EIO;
  if (bytesWritten < 0) return errno;
  
  return 0;
}

int ext_bio_writeDouble(ext_binaryIO* bio, double d)
{
  if (bio == NULL) return EFAULT;
  
  ssize_t bytesWritten = write(bio->fileDescriptor, &d, sizeof(double));
  if (bytesWritten == 0) return EIO;
  if (bytesWritten < 0) return errno;
  
  return 0;
} 

int ext_bio_writeDoubles(ext_binaryIO* bio, const double* d, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  int errorCode = ext_bio_writeSizeType(bio, length);
  if (errorCode != 0) return errorCode;
  
  return ext_bio_writeNDoubles(bio, d, length);
}

int ext_bio_writeNDoubles(ext_binaryIO* bio, const double* d, size_t length)
{
  size_t totalItemsWritten = 0;
  while (totalItemsWritten < length) {
    size_t itemsWritten = fillBufferFromDoubles(bio, d + totalItemsWritten, length - totalItemsWritten);
    
    size_t totalBytesWritten = 0;
    size_t bytesToWrite = itemsWritten * sizeof(double);
    
    while (totalBytesWritten < bytesToWrite) {
      ssize_t bytesWritten = write(bio->fileDescriptor, (char*) bio->buffer, bytesToWrite - totalBytesWritten);
      if (bytesWritten == 0) return EIO;
      if (bytesWritten < 0) return errno;
      totalBytesWritten += (size_t) bytesWritten;
    }
    
    totalItemsWritten += itemsWritten;
  }
  
  return 0;
}

// these return NULL if there is an error and set errno appropriately
int ext_bio_readChar(ext_binaryIO* bio, char *c)
{
  if (bio == NULL) return EFAULT;
  
  ssize_t bytesRead = read(bio->fileDescriptor, c, sizeof(char));
  
  if (bytesRead == 0) return EIO;
  if (bytesRead < 0) return errno;
  
  return 0;
}

char* ext_bio_readChars(ext_binaryIO* bio, size_t* length)
{
  if (bio == NULL) { errno = EFAULT; return NULL; }
  
  int errorCode = ext_bio_readSizeType(bio, length);
  if (errorCode != 0) { errno = errorCode; return NULL; }
  
  char* c = (char*) malloc(*length * sizeof(char));
  if (c == NULL) { errno = ENOMEM; return NULL; }
  
  errorCode = ext_bio_readNChars(bio, c, *length);
  if (errorCode != 0) {
    free(c);
    errno = errorCode;
    return NULL;
  }
  
  return c;
}

int ext_bio_readNChars(ext_binaryIO* bio, char* c, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  size_t totalBytesRead = 0;
  while (totalBytesRead < length) {
    ssize_t bytesRead = read(bio->fileDescriptor, c + totalBytesRead, length - totalBytesRead);
    if (bytesRead == 0) return EIO;
    if (bytesRead < 0) return errno;
    totalBytesRead += (size_t) bytesRead;
  }
  
  return 0;
}


int ext_bio_readSizeType(ext_binaryIO* bio, size_t* s)
{
  if (bio == NULL) return EFAULT;
  
  uint64_t u;
  
  ssize_t bytesRead = read(bio->fileDescriptor, &u, sizeof(uint64_t));
  if (bytesRead == 0) return EIO;
  if (bytesRead < 0) return errno;
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor8ByteWord((char*) &u);
#endif
  
  if (u > (uint64_t) SIZE_MAX) return EOVERFLOW;
  
  *s = (size_t) u;
  
  return 0;
}

size_t* ext_bio_readSizeTypes(ext_binaryIO* bio, size_t* length)
{
  if (bio == NULL) { errno = EFAULT; return NULL; }
  
  int errorCode = ext_bio_readSizeType(bio, length);
  if (errorCode != 0) { errno = errorCode; return NULL; }
  
  size_t* s = (size_t*) malloc(*length * sizeof(size_t));
  if (s == NULL) { errno = ENOMEM; return NULL; }
  
  errorCode = ext_bio_readNSizeTypes(bio, s, *length);
  if (errorCode != 0) {
    free(s);
    errno = errorCode;
    return NULL;
  }
  
  return s;
}

int ext_bio_readNSizeTypes(ext_binaryIO* bio, size_t* s, size_t length)
{
  if (bio == NULL) return errno = EFAULT;
  
  size_t totalItemsRead = 0;
  while (totalItemsRead < length) {
    size_t itemsToRead = (bio->bufferLength / sizeof(uint64_t)) < (length - totalItemsRead) ? bio->bufferLength / sizeof(uint64_t) : length - totalItemsRead;
    
    size_t totalBytesRead = 0;
    size_t bytesToRead = itemsToRead * sizeof(uint64_t);
    
    while (totalBytesRead < bytesToRead) {
      ssize_t bytesRead = read(bio->fileDescriptor, (char*) bio->buffer, bytesToRead - totalBytesRead);
      if (bytesRead == 0) return EIO;
      if (bytesRead < 0) return errno;
      totalBytesRead += (size_t) bytesRead;
    }
    
    size_t itemsRead = fillSizeTypesFromBuffer(bio, s + totalItemsRead, itemsToRead);
    if (itemsRead == 0) return errno;
    
    totalItemsRead += itemsRead;
  }
  
  return 0;
}

int ext_bio_readUnsigned32BitInteger(ext_binaryIO* bio, uint32_t* u)
{
  if (bio == NULL) return EFAULT;
  
  ssize_t bytesRead = read(bio->fileDescriptor, u, sizeof(uint32_t));
  if (bytesRead == 0) return EIO;
  if (bytesRead < 0) return errno;
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor4ByteWord((char*) u);
#endif
  
  return 0;
}

uint32_t* ext_bio_readUnsigned32BitIntegers(ext_binaryIO* bio, size_t* length)
{
  if (bio == NULL) { errno = EFAULT; return NULL; }
  
  int errorCode = ext_bio_readSizeType(bio, length);
  if (errorCode != 0) { errno = errorCode; return NULL; }
  
  uint32_t* u = (uint32_t*) malloc(*length * sizeof(uint32_t));
  if (u == NULL) { errno = ENOMEM; return NULL; }
  
  errorCode = ext_bio_readNUnsigned32BitIntegers(bio, u, *length);
  if (errorCode != 0) {
    free(u);
    errno = errorCode;
    return NULL;
  }
  
  return 0;
}

int ext_bio_readNUnsigned32BitIntegers(ext_binaryIO* bio, uint32_t* u, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  size_t totalItemsRead = 0;
  while (totalItemsRead < length) {
    size_t itemsToRead = (bio->bufferLength / sizeof(uint32_t)) < (length - totalItemsRead) ? bio->bufferLength / sizeof(uint32_t) : length - totalItemsRead;
    
    size_t totalBytesRead = 0;
    size_t bytesToRead = itemsToRead * sizeof(uint32_t);
    
    while (totalBytesRead < bytesToRead) {
      ssize_t bytesRead = read(bio->fileDescriptor, (char*) bio->buffer, bytesToRead - totalBytesRead);
      if (bytesRead == 0) return EIO;
      if (bytesRead < 0) return errno;
      totalBytesRead += (size_t) bytesRead;
    }
    
    totalItemsRead += fillUnsigned32BitIntegersFromBuffer(bio, u + totalItemsRead, itemsToRead);
  }
  
  return 0;
}

int ext_bio_readUnsigned64BitInteger(ext_binaryIO* bio, uint64_t* u)
{
  if (bio == NULL) return EFAULT;
  
  ssize_t bytesRead = read(bio->fileDescriptor, u, sizeof(uint64_t));
  if (bytesRead == 0) return EIO;
  if (bytesRead < 0) return errno;
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor8ByteWord((char*) u);
#endif
  
  return 0;
}

int ext_bio_readDouble(ext_binaryIO* bio, double* d)
{
  if (bio == NULL) return EFAULT;
  
  ssize_t bytesRead = read(bio->fileDescriptor, d, sizeof(double));
  if (bytesRead == 0) return EIO;
  if (bytesRead < 0) return errno;
    
  return 0;
}

double* ext_bio_readDoubles(ext_binaryIO* bio, size_t* length)
{
  if (bio == NULL) { errno = EFAULT; return NULL; }
  
  int errorCode = ext_bio_readSizeType(bio, length);
  if (errorCode != 0) { errno = errorCode; return NULL; }
  
  double* d = (double*) malloc(*length * sizeof(double));
  if (d == NULL) { errno = ENOMEM; return NULL; }
  
  errorCode = ext_bio_readNDoubles(bio, d, *length);
  if (errorCode != 0) {
    free(d);
    errno = errorCode;
    return NULL;
  }
  
  return 0;
}

int ext_bio_readNDoubles(ext_binaryIO* bio, double* d, size_t length)
{
  if (bio == NULL) return EFAULT;
  
  size_t totalItemsRead = 0;
  while (totalItemsRead < length) {
    size_t itemsToRead = (bio->bufferLength / sizeof(double)) < (length - totalItemsRead) ? bio->bufferLength / sizeof(double) : length - totalItemsRead;
    
    size_t totalBytesRead = 0;
    size_t bytesToRead = itemsToRead * sizeof(double);
    
    while (totalBytesRead < bytesToRead) {
      ssize_t bytesRead = read(bio->fileDescriptor, (char*) bio->buffer, bytesToRead - totalBytesRead);
      if (bytesRead == 0) return EIO;
      if (bytesRead < 0) return errno;
      totalBytesRead += (size_t) bytesRead;
    }
    
    totalItemsRead += fillDoublesFromBuffer(bio, d + totalItemsRead, itemsToRead);
  }
  
  return 0;
}

// upcast to 64_bit ints
static size_t fillBufferFromSizeTypes(ext_binaryIO* restrict bio, const size_t* restrict v, size_t length)
{
  if (length == 0) return 0;
  
  // purposefully aligned this stupid pointer
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  uint64_t* restrict buffer = (uint64_t* restrict) bio->buffer;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  size_t bufferLength = bio->bufferLength / sizeof(uint64_t);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
  if (sizeof(size_t) == sizeof(uint64_t)) {
    memcpy(buffer, v, fillLength * sizeof(uint64_t));
#ifndef WORDS_BIGENDIAN
   swapEndiannessFor8ByteWords((char*) buffer, fillLength);
#endif
    return fillLength;
  }

  size_t lengthMod5 = fillLength % 5;
  size_t i = 0;
  for ( ; i < lengthMod5; ++i) buffer[i] = (uint64_t) v[i];
  
  for ( ; i < fillLength; i += 5) {
    buffer[    i] = (uint64_t) v[i];
    buffer[i + 1] = (uint64_t) v[i + 1];
    buffer[i + 2] = (uint64_t) v[i + 2];
    buffer[i + 3] = (uint64_t) v[i + 3];
    buffer[i + 4] = (uint64_t) v[i + 4];
  }
  
#ifndef WORDS_BIGENDIAN
   swapEndiannessFor8ByteWords((char*) buffer, fillLength);
#endif
  
  return fillLength;
}

// upcast to 64_bit ints
static size_t fillSizeTypesFromBuffer(ext_binaryIO* restrict bio, size_t* restrict v, size_t length)
{
  if (length == 0) return 0;
  
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  uint64_t* restrict buffer = (uint64_t* restrict) bio->buffer;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  size_t bufferLength = bio->bufferLength / sizeof(uint64_t);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor8ByteWords((char*) buffer, fillLength);
#endif
  
  if (sizeof(size_t) == sizeof(uint64_t)) {
    memcpy(v, buffer, fillLength * sizeof(uint64_t));
    
    return fillLength;
  }
  
  size_t i = 0;
  for ( ; i < fillLength; ++i) {
    if (buffer[i] > (uint64_t) SIZE_MAX) {
      errno = EOVERFLOW;
      return 0;
    }
  } 
  size_t lengthMod5 = fillLength % 5;
  for ( ; i < lengthMod5; ++i) v[i] = (size_t) buffer[i];
  
  for ( ; i < fillLength; i += 5) {
    v[    i] = (size_t) buffer[i];
    v[i + 1] = (size_t) buffer[i + 1];
    v[i + 2] = (size_t) buffer[i + 2];
    v[i + 3] = (size_t) buffer[i + 3];
    v[i + 4] = (size_t) buffer[i + 4];
  }
  
  return fillLength;
}

static size_t fillBufferFromDoubles(ext_binaryIO* restrict bio, const double* restrict d, size_t length)
{
  if (length == 0) return 0;
  
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  double* restrict buffer = (double* restrict) bio->buffer;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  size_t bufferLength = bio->bufferLength / sizeof(double);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
  memcpy(buffer, d, fillLength * sizeof(double));
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor8ByteWords((char*) buffer, fillLength);
#endif
   
  return fillLength;
}

static size_t fillDoublesFromBuffer(ext_binaryIO* restrict bio, double* restrict d, size_t length)
{
  if (length == 0) return 0;
  
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  double* restrict buffer = (double* restrict) bio->buffer;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  size_t bufferLength = bio->bufferLength / sizeof(double);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor8ByteWords((char*) buffer, fillLength);
#endif
  
  memcpy(d, buffer, fillLength * sizeof(double));
   
  return fillLength;
}

static size_t fillBufferFromUnsigned32BitIntegers(ext_binaryIO* restrict bio, const uint32_t* restrict v, size_t length)
{
  if (length == 0) return 0;
  
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  uint32_t* restrict buffer = (uint32_t* restrict) bio->buffer;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  size_t bufferLength = bio->bufferLength / sizeof(uint32_t);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
  memcpy(buffer, v, fillLength * sizeof(uint32_t));
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor4ByteWords((char*) buffer, fillLength);
#endif
   
  return fillLength;
}

static size_t fillUnsigned32BitIntegersFromBuffer(ext_binaryIO* restrict bio, uint32_t* restrict v, size_t length)
{
  if (length == 0) return 0;
  
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  uint32_t* restrict buffer = (uint32_t* restrict) bio->buffer;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  size_t bufferLength = bio->bufferLength / sizeof(uint32_t);
  size_t fillLength = (length < bufferLength ? length : bufferLength);
  
#ifndef WORDS_BIGENDIAN
  swapEndiannessFor4ByteWords((char*) buffer, fillLength);
#endif
  
  memcpy(v, buffer, fillLength * sizeof(uint32_t));
   
  return fillLength;
}

#ifndef WORDS_BIGENDIAN
static void swapEndiannessFor4ByteWords(char* c, size_t length)
{
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  uint32_t* u = (uint32_t*) c;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  
  size_t lengthMod5 = length % 5;
  size_t i = 0;
  for ( ; i < lengthMod5; ++i) swapEndiannessFor4ByteWord((char*) (u + i));
  
  for ( ; i < length; i += 5) {
    swapEndiannessFor4ByteWord((char*) (u + i));
    swapEndiannessFor4ByteWord((char*) (u + i + 1));
    swapEndiannessFor4ByteWord((char*) (u + i + 2));
    swapEndiannessFor4ByteWord((char*) (u + i + 3));
    swapEndiannessFor4ByteWord((char*) (u + i + 4));
  }
}

static void swapEndiannessFor8ByteWords(char* c, size_t length)
{
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-align"
#endif
  uint64_t* u = (uint64_t*) c;
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#  pragma GCC diagnostic pop
#endif
  
  size_t lengthMod5 = length % 5;
  size_t i = 0;
  for ( ; i < lengthMod5; ++i) swapEndiannessFor8ByteWord((char*) (u + i));
  
  for ( ; i < length; i += 5) {
    swapEndiannessFor8ByteWord((char*) (u + i));
    swapEndiannessFor8ByteWord((char*) (u + i + 1));
    swapEndiannessFor8ByteWord((char*) (u + i + 2));
    swapEndiannessFor8ByteWord((char*) (u + i + 3));
    swapEndiannessFor8ByteWord((char*) (u + i + 4));
  }
}

#endif
