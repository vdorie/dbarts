#include <external/thread.h>
#include "pthread.c"

#include "config.h"

#include <stdlib.h>
#include <errno.h>
#include <stdbool.h>

// clock_gettime + CLOCK_REALTIME are in time.h, gettimeofday is in sys/time.h; plain time() is in time.h too
#if (!defined(HAVE_CLOCK_GETTIME) || !defined(CLOCK_REALTIME)) && defined(HAVE_GETTIMEOFDAY)
#  include <sys/time.h>
#endif
#include <time.h>


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

struct ThreadData;

static int initializeManager(ext_mt_manager_t manager, size_t numThreads);
static int initializeThreadData(ext_mt_manager_t manager, struct ThreadData* data, size_t index);
static int destroyThreadData(struct ThreadData* data);

static void* threadLoop(void* _data);

typedef struct ThreadData {
  ext_mt_manager_t manager;
  Condition taskAvailable;
  size_t index;
  
  ext_mt_taskFunction_t task;
  void* taskData;
} ThreadData;

typedef struct {
  size_t* elements;
  size_t queueSize;
  size_t pushIndex;
  size_t popIndex;
} IndexArrayQueue;

UNUSED static IndexArrayQueue* createIndexArrayQueue(size_t queueSize); // returns NULL if an error occurs
UNUSED static void destroyIndexArrayQueue(IndexArrayQueue* queue);
static int initializeIndexArrayQueue(IndexArrayQueue* queue, size_t queueSize);
static void invalidateIndexArrayQueue(IndexArrayQueue* queue);

static size_t IAQ_INVALID = ((size_t) -1);
static int push(IndexArrayQueue* queue, size_t element); // can return ENOBUFS if full
static size_t pop(IndexArrayQueue* queue); // returns IAQ_INVALID if buffer is empty
static size_t getNumElementsInQueue(const IndexArrayQueue* queue);

typedef struct _ext_mt_manager_t {
  Thread* threads;
  ThreadData* threadData;
  
  IndexArrayQueue threadQueue;
  
  size_t numThreads;
  size_t numThreadsActive;  // active means created, but possibly waiting
  size_t numThreadsRunning; // running means actually computing some stuff
  
  bool threadsShouldExit;
  Mutex mutex;
  Condition threadIsActive;
  Condition taskDone;
} _ext_mt_manager_t;


int ext_mt_create(ext_mt_manager_t* managerPtr, size_t numThreads)
{
  *managerPtr = (ext_mt_manager_t) malloc(sizeof(_ext_mt_manager_t));
  if (*managerPtr == NULL) return ENOMEM;
  
  ext_mt_manager_t manager = *managerPtr;
  
  int result = initializeManager(manager, numThreads);
  if (result != 0) {
    *managerPtr = NULL;
    return result;
  }
  
  for (size_t i = 0; i < numThreads; ++i) {
    result = initializeThreadData(manager, &manager->threadData[i], i);
    
    if (result != 0) { manager->numThreads = i; break; }
    
    result = createThread(&manager->threads[i], &threadLoop, &manager->threadData[i]);
    
    if (result != 0) {
      destroyThreadData(&manager->threadData[i]);
      manager->numThreads = i;
      break;
    }
    
    push(&manager->threadQueue, i);
  }
  
  // wait for all threads to check in
  lockMutex(manager->mutex);
  while (manager->numThreadsActive < manager->numThreads)
    waitOnCondition(manager->threadIsActive, manager->mutex);
  unlockMutex(manager->mutex);
  
  if (result != 0) {
    ext_mt_destroy(manager);
    *managerPtr = NULL;
  }
  
  return result;
}

/* size_t ext_mt_getThreadId(const ext_mt_manager_t manager)
{
  Thread nativeThreadId = pthread_self();
  size_t i;
  for (i = 0; i < manager->numThreads; ++i) if (nativeThreadId == manager->threads[i]) break;
  
  return i;
} */

size_t ext_mt_getNumThreads(const ext_mt_manager_t manager)
{
  if (manager == NULL) return 1;
  return manager->numThreads;
}


void ext_mt_getNumThreadsForJob(const ext_mt_manager_t restrict threadManager, size_t numElements, size_t minNumElementsPerThread,
                                size_t* restrict numThreadsPtr, size_t* restrict numElementsPerThreadPtr, size_t* restrict offByOneIndexPtr)
{
  size_t numThreadsManaged = 0;
  if (numElements < 2 * minNumElementsPerThread || threadManager == NULL ||
      (numThreadsManaged = ext_mt_getNumThreads(threadManager)) <= 1) {
    if (numThreadsPtr != NULL) *numThreadsPtr = 1;
    *numElementsPerThreadPtr = numElements;
    *offByOneIndexPtr = 1;
    return;
  }
  
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

int ext_mt_runTasks(ext_mt_manager_t restrict manager, ext_mt_taskFunction_t function,
                    void** restrict data, size_t numTasks)
{
  if (manager->threads == NULL || manager->threadData == NULL ||
      manager->numThreadsActive == 0) return EINVAL;
  
  int result = 0;
  
  ThreadData* threadData = manager->threadData;
  
  lockMutex(manager->mutex);
  
  for (size_t i = 0; i < numTasks; ++i) {
    while (getNumElementsInQueue(&manager->threadQueue) == 0) waitOnCondition(manager->taskDone, manager->mutex);
    
    size_t j = pop(&manager->threadQueue);
    
    threadData[j].task = function;
    threadData[j].taskData = (data == NULL ? NULL : data[i]);
    manager->numThreadsRunning++;
      
    signalCondition(threadData[j].taskAvailable);
  }
  
  while (manager->numThreadsRunning > 0) waitOnCondition(manager->taskDone, manager->mutex);
  
  unlockMutex(manager->mutex);
  
  return result;
}

static void inline getTime(struct timespec* ts)
{
#if defined(HAVE_CLOCK_GETTIME) && defined(CLOCK_REALTIME)
  clock_gettime(CLOCK_REALTIME, ts);
#elif defined(HAVE_GETTIMEOFDAY)
  struct timeval tv;
  gettimeofday(&tv, NULL);
  ts->tv_sec = tv.tv_sec;
  ts->tv_nsec = 1000 * tv.tv_usec;
#else
  ts->tv_sec = time(NULL);
  ts->tv_nsec = 0;
#endif
}

int ext_mt_runTasksWithInfo(ext_mt_manager_t restrict manager, ext_mt_taskFunction_t function,
                            void** restrict data, size_t numTasks, time_t sleepSeconds, ext_mt_infoFunction_t info)
{
  if (manager->threads == NULL || manager->threadData == NULL ||
      manager->numThreadsActive == 0) return EINVAL;
  
  int result = 0;
  
  ThreadData* threadData = manager->threadData;
  
  struct timespec wakeTime;
  
  lockMutex(manager->mutex);
  
  getTime(&wakeTime);
    
  wakeTime.tv_sec += sleepSeconds;
  for (size_t i = 0; i < numTasks; ++i) {
    // while (getNumElementsInQueue(&manager->threadQueue) == 0) waitOnCondition(manager->taskDone, manager->mutex);
    
    while (getNumElementsInQueue(&manager->threadQueue) == 0) {
      int waitStatus = waitOnConditionForTime(manager->taskDone, manager->mutex, wakeTime);
      if (waitStatus == ETIMEDOUT) {
        if (info != NULL) info(data, manager->numThreads);
        
        getTime(&wakeTime);
        wakeTime.tv_sec += sleepSeconds;
      }
    }
    
    size_t j = pop(&manager->threadQueue);
    
    threadData[j].task = function;
    threadData[j].taskData = (data == NULL ? NULL : data[i]);
    manager->numThreadsRunning++;
    
    signalCondition(threadData[j].taskAvailable);
  }
  
  while (manager->numThreadsRunning > 0) {
    int waitStatus = waitOnConditionForTime(manager->taskDone, manager->mutex, wakeTime);
    if (waitStatus == ETIMEDOUT) {
      if (info != NULL) info(data, manager->numThreads);
      
      getTime(&wakeTime);
      wakeTime.tv_sec += sleepSeconds;
    }
  }
  
  unlockMutex(manager->mutex);
  
  return result;
}


static void* threadLoop(void* _data)
{
  ThreadData* data = (ThreadData*) _data;
  ext_mt_manager_t manager = data->manager;
  
  lockMutex(manager->mutex);
  manager->numThreadsActive++;
  signalCondition(manager->threadIsActive);
  
  
  while (true) {
    while (data->task == NULL && manager->threadsShouldExit == false) waitOnCondition(data->taskAvailable, manager->mutex);
    if (manager->threadsShouldExit == true) break;
    
    unlockMutex(manager->mutex);
    
    data->task(data->taskData);
    
    lockMutex(manager->mutex);
    manager->numThreadsRunning--;
    data->task = NULL;
    data->taskData = NULL;
    signalCondition(manager->taskDone);
    push(&manager->threadQueue, data->index);
  }
  
  manager->numThreadsActive--;
  
  unlockMutex(manager->mutex);
  
  return NULL;
}

int ext_mt_destroy(ext_mt_manager_t manager)
{
  if (manager == NULL) return 0;
    
  int result = 0;
  
  if (manager->threads != NULL && manager->threadData != NULL &&
      manager->numThreadsActive > 0 && manager->numThreads > 0) {
    ThreadData* threadData = manager->threadData;
    
    lockMutex(manager->mutex);
    manager->threadsShouldExit = true;
    for (size_t i = 0; i < manager->numThreads; ++i)
      signalCondition(threadData[i].taskAvailable);
    
    unlockMutex(manager->mutex);
    
    for (size_t i = 0; i < manager->numThreads; ++i)
      result |= joinThread(manager->threads[i]);
  }
  
  invalidateIndexArrayQueue(&manager->threadQueue);
  
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  
  if (manager->threadData != NULL) {
    for (size_t i = 0; i < manager->numThreads; ++i) result |= destroyThreadData(&manager->threadData[i]);
    
    free(manager->threadData);
    manager->threadData = NULL;
  }
  
  result |= destroyMutex(manager->mutex);
  result |= destroyCondition(manager->threadIsActive);
  result |= destroyCondition(manager->taskDone);
  
  free(manager);
  
  return result;
}

static int initializeManager(ext_mt_manager_t manager, size_t numThreads)
{
  int result;
  
  manager->numThreadsRunning = 0;
  manager->numThreadsActive = 0;
  manager->numThreads = numThreads;
  manager->threads = NULL;
  manager->threadData = NULL;
  manager->threadsShouldExit = false;
  bool mutexInitialized = false;
  bool threadIsActiveInitialized = false;
  bool taskDoneInitialized = false;
  
  manager->threads = (Thread*) malloc(numThreads * sizeof(Thread));
  if (manager->threads == NULL) { result = ENOMEM; goto ext_mt_initialization_failed; }
  
  manager->threadData = (ThreadData*) malloc(numThreads * sizeof(ThreadData));
  if (manager->threadData == NULL) { result = ENOMEM; goto ext_mt_initialization_failed; }
    
  result = initializeIndexArrayQueue(&manager->threadQueue, numThreads);
  if (result != 0) goto ext_mt_initialization_failed;
  
  result = initializeMutex(manager->mutex);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) mutexInitialized = true;
    goto ext_mt_initialization_failed;
  }
  
  result = initializeCondition(manager->threadIsActive);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) threadIsActiveInitialized = true;
    goto ext_mt_initialization_failed;
  }
  
  result = initializeCondition(manager->taskDone);
  if (result != 0) {
    taskDoneInitialized = true;
    goto ext_mt_initialization_failed;
  }
  
  return result;
  
ext_mt_initialization_failed:
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  if (manager->threadData != NULL) { free(manager->threadData); manager->threadData = NULL; }
  
  invalidateIndexArrayQueue(&manager->threadQueue);
  
  if (mutexInitialized) destroyMutex(manager->mutex);
  if (threadIsActiveInitialized) destroyCondition(manager->threadIsActive);
  if (taskDoneInitialized) destroyCondition(manager->taskDone);
  
  free(manager);
  
  return result;
}

static int initializeThreadData(ext_mt_manager_t manager, ThreadData* data, size_t index)
{
  data->manager = manager;
  data->index = index;
  
  data->task = NULL;
  data->taskData = NULL;
  
  int result = initializeCondition(data->taskAvailable);
  
  if (result != 0 && result != EBUSY && result != EINVAL)
    destroyCondition(data->taskAvailable);
  
  return result;
}

static int destroyThreadData(ThreadData* data)
{
  return destroyCondition(data->taskAvailable);
}

static IndexArrayQueue* createIndexArrayQueue(size_t queueSize)
{
  IndexArrayQueue* result = (IndexArrayQueue*) malloc(sizeof(IndexArrayQueue));
  if (result == NULL) return result;
  
  if (initializeIndexArrayQueue(result, queueSize) != 0) {
    free(result);
    return NULL;
  }
  
  return result;
}

static int initializeIndexArrayQueue(IndexArrayQueue* queue, size_t queueSize)
{
  queue->elements = (size_t*) malloc(queueSize * sizeof(size_t));
  if (queue->elements == NULL) return ENOMEM;

  for (size_t i = 0; i < queueSize; ++i) queue->elements[i] = IAQ_INVALID;
  
  queue->queueSize = queueSize;
  queue->pushIndex = 0;
  queue->popIndex = 0;
  
  return 0;
}

static void destroyIndexArrayQueue(IndexArrayQueue* queue)
{
  invalidateIndexArrayQueue(queue);
  if (queue != NULL) free(queue);
}

static void invalidateIndexArrayQueue(IndexArrayQueue* queue)
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

static int push(IndexArrayQueue* queue, size_t element)
{
  if (queue->pushIndex == queue->popIndex && queue->elements[queue->pushIndex] != IAQ_INVALID) return ENOBUFS;
                                                             
  queue->elements[queue->pushIndex++] = element;
  if (queue->pushIndex == queue->queueSize) queue->pushIndex = 0;
  
  return 0;
}

static size_t pop(IndexArrayQueue* queue)
{
  if (queue->popIndex == queue->pushIndex && queue->elements[queue->popIndex] == IAQ_INVALID) return IAQ_INVALID;
  
  size_t result = queue->elements[queue->popIndex];
  queue->elements[queue->popIndex++] = IAQ_INVALID;
  if (queue->popIndex == queue->queueSize) queue->popIndex = 0;
  
  return result;
}

static size_t getNumElementsInQueue(const IndexArrayQueue* queue)
{
  if (queue->popIndex == queue->pushIndex) {
    if (queue->elements[queue->popIndex] == IAQ_INVALID) return 0;
    return queue->queueSize;
  }
  
  size_t pushIndex = queue->pushIndex < queue->popIndex ? queue->pushIndex + queue->queueSize  : queue->pushIndex;
  
  return pushIndex - queue->popIndex;
}
