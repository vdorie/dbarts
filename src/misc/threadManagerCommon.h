#ifndef MISC_THREAD_MANAGER_COMMON_H
#define MISC_THREAD_MANAGER_COMMON_H

// Internal helper for the thread-manager translation unit (thread.c). This is
// a private header - it is NOT part of the public misc/thread.h API.
//
// misc_partitionThreadJob is the integer split arithmetic the manager runs to
// hand out elements across its threads. (Its ready queue is the ring buffer in
// indexArrayQueue.h, likewise private.)

#include <stddef.h> // size_t

// Split numElements across numThreadsManaged threads/pieces, assigning
// numElementsPerThread to the first offByOneIndex threads and one fewer to the
// rest. Callers handle the trivial (single-thread) short-circuit themselves;
// this runs once numThreadsManaged (> 1) has been established.
static inline void misc_partitionThreadJob(size_t numElements, size_t minNumElementsPerThread, size_t numThreadsManaged,
                                           size_t* restrict numThreadsPtr, size_t* restrict numElementsPerThreadPtr,
                                           size_t* restrict offByOneIndexPtr)
{
  size_t numThreads = minNumElementsPerThread == 0 ? numElements : numElements / minNumElementsPerThread;
  if (numThreads > numThreadsManaged) numThreads = numThreadsManaged;

  size_t numElementsPerThread = numElements / numThreads;
  size_t offByOneIndex = numElements % numThreads;
  if (offByOneIndex != 0) numElementsPerThread += 1;
  else offByOneIndex = numThreads;

  if (numThreadsPtr != NULL) *numThreadsPtr = numThreads;
  *numElementsPerThreadPtr = numElementsPerThread;
  *offByOneIndexPtr = offByOneIndex;
}

#endif // MISC_THREAD_MANAGER_COMMON_H
