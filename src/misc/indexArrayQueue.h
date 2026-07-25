#ifndef MISC_INDEX_ARRAY_QUEUE_H
#define MISC_INDEX_ARRAY_QUEUE_H

// A fixed-capacity ring buffer of size_t indices, shared verbatim by the plain
// (thread.c) and blocking (blockingThreadManager.c) thread managers. Private
// header - NOT part of the public misc/thread.h API. Empty slots hold
// IAQ_INVALID; push returns ENOBUFS when full, pop returns IAQ_INVALID when
// empty.
//
// Deliberately NOT included by hierarchicalThreadManager.c: that unit has its
// own push/pop over a linked ThreadStack and does not use this queue.

#include <stddef.h> // size_t
#include <errno.h>  // ENOMEM, ENOBUFS
#include <stdlib.h> // malloc, free

typedef struct {
  size_t* elements;
  size_t queueSize;
  size_t pushIndex;
  size_t popIndex;
} IndexArrayQueue;

static size_t IAQ_INVALID = ((size_t) -1);

static inline int initializeIndexArrayQueue(IndexArrayQueue* queue, size_t queueSize)
{
  queue->elements = (size_t*) malloc(queueSize * sizeof(size_t));
  if (queue->elements == NULL) return ENOMEM;

  for (size_t i = 0; i < queueSize; ++i) queue->elements[i] = IAQ_INVALID;

  queue->queueSize = queueSize;
  queue->pushIndex = 0;
  queue->popIndex = 0;

  return 0;
}

static inline void invalidateIndexArrayQueue(IndexArrayQueue* queue)
{
  if (queue == NULL || queue->elements == NULL) return;

  free(queue->elements);
  queue->elements = NULL;
}

#ifndef ENOBUFS
#ifdef __WIN32
#include <error.h>
#define ENOBUFS ERROR_INSUFFICIENT_BUFFER
#else
#define ENOBUFS 105
#endif
#endif

// can return ENOBUFS if full
static inline int push(IndexArrayQueue* queue, size_t element)
{
  if (queue->pushIndex == queue->popIndex && queue->elements[queue->pushIndex] != IAQ_INVALID) return ENOBUFS;

  queue->elements[queue->pushIndex++] = element;
  if (queue->pushIndex == queue->queueSize) queue->pushIndex = 0;

  return 0;
}

// returns IAQ_INVALID if buffer is empty
static inline size_t pop(IndexArrayQueue* queue)
{
  if (queue->popIndex == queue->pushIndex && queue->elements[queue->popIndex] == IAQ_INVALID) return IAQ_INVALID;

  size_t result = queue->elements[queue->popIndex];
  queue->elements[queue->popIndex++] = IAQ_INVALID;
  if (queue->popIndex == queue->queueSize) queue->popIndex = 0;

  return result;
}

static inline size_t getNumElementsInQueue(const IndexArrayQueue* queue)
{
  if (queue->popIndex == queue->pushIndex) {
    if (queue->elements[queue->popIndex] == IAQ_INVALID) return 0;
    return queue->queueSize;
  }

  size_t pushIndex = queue->pushIndex < queue->popIndex ? queue->pushIndex + queue->queueSize  : queue->pushIndex;

  return pushIndex - queue->popIndex;
}

#endif // MISC_INDEX_ARRAY_QUEUE_H
