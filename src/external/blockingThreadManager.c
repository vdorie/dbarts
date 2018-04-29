#include <external/thread.h>
#include "pthread.c"

#include "config.h"

#include <stdlib.h>
#include <errno.h>
#include <stdbool.h>

// clock_gettime + CLOCK_REALTIME are in time.h, gettimeofday is in sys/time.h; plain time() is in time.h too
// time.h imported from <external/thread.h>
#if (!defined(HAVE_CLOCK_GETTIME) || !defined(CLOCK_REALTIME)) && defined(HAVE_GETTIMEOFDAY)
#  include <sys/time.h>
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

struct ThreadData;

static int initializeManager(ext_btm_manager_t manager, size_t numThreads);
static int initializeThreadData(ext_btm_manager_t manager, struct ThreadData* data, size_t index);
static int destroyThreadData(struct ThreadData* data);

static void* threadLoop(void* _data);

typedef struct {
  void** elements;
  size_t queueSize;
  size_t pushIndex;
  size_t popIndex;
} ArrayQueue;

typedef struct {
  size_t* elements;
  size_t queueSize;
  size_t pushIndex;
  size_t popIndex;
} IndexArrayQueue;


static size_t IAQ_INVALID = ((size_t) -1);
static int initializeIndexArrayQueue(IndexArrayQueue* queue, size_t queueSize);
static void invalidateIndexArrayQueue(IndexArrayQueue* queue);

static int push(IndexArrayQueue* queue, size_t element); // can return ENOBUFS if full
static size_t pop(IndexArrayQueue* queue); // returns IAQ_INVALID if buffer is empty
static size_t getNumElementsInQueue(const IndexArrayQueue* queue);

typedef struct ThreadData {
  ext_btm_manager_t manager;
  Condition taskAvailable;
  Condition parentTaskComplete;
  size_t id;
  
  ext_btm_taskFunction_t task;
  void* taskData;
  
  ext_btm_taskFunction_t parentTask;
  void* parentTaskData;
  bool parentIsFinished;
} ThreadData;

typedef struct _ext_btm_manager_t {
  Thread* threads;
  ThreadData* threadData;
  
  IndexArrayQueue threadQueue;
  IndexArrayQueue parentTaskQueue;
  
  size_t numThreads;
  size_t numThreadsActive;  // active means created, but possibly waiting
  size_t numThreadsRunning; // running means actually computing some stuff
  
  bool threadsShouldExit;
  Mutex mutex;
  Condition threadIsActive;
  Condition threadIsWaiting;
  
} _ext_btm_manager_t;

#include <stdio.h>

int ext_btm_create(ext_btm_manager_t* managerPtr, size_t numThreads)
{
  *managerPtr = (ext_btm_manager_t) malloc(sizeof(_ext_btm_manager_t));
  if (*managerPtr == NULL) return ENOMEM;
  
  ext_btm_manager_t manager = *managerPtr;
  
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
    ext_btm_destroy(manager);
    *managerPtr = NULL;
  }
  
  return result;
}

size_t ext_btm_getThreadId(const ext_btm_manager_t manager)
{
  Thread nativeThreadId = pthread_self();
  size_t i;
  for (i = 0; i < manager->numThreads; ++i) if (nativeThreadId == manager->threads[i]) break;
  
  return i;
}

size_t ext_btm_getNumThreads(const ext_btm_manager_t manager)
{
  if (manager == NULL) return 1;
  return manager->numThreads;
}


void ext_btm_getNumThreadsForJob(const ext_btm_manager_t restrict threadManager, size_t numElements, size_t minNumElementsPerThread,
                                size_t* restrict numThreadsPtr, size_t* restrict numElementsPerThreadPtr, size_t* restrict offByOneIndexPtr)
{
  size_t numThreadsManaged = 0;
  if (numElements < 2 * minNumElementsPerThread || threadManager == NULL ||
      (numThreadsManaged = ext_btm_getNumThreads(threadManager)) <= 1) {
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

int ext_btm_runTasks(ext_btm_manager_t restrict manager, ext_btm_taskFunction_t function,
                     void** restrict data, size_t numTasks)
{
  if (manager->threads == NULL || manager->threadData == NULL ||
      manager->numThreadsActive == 0) return EINVAL;
  
  int result = 0;
  
  ThreadData* threadData = manager->threadData;
  
  lockMutex(manager->mutex);
  
  for (size_t i = 0; i < numTasks; /* */ ) {
    while (getNumElementsInQueue(&manager->threadQueue) == 0 && getNumElementsInQueue(&manager->parentTaskQueue) == 0)
      waitOnCondition(manager->threadIsWaiting, manager->mutex);
    
    while (getNumElementsInQueue(&manager->parentTaskQueue) != 0) {
      size_t j = pop(&manager->parentTaskQueue);
      
      threadData[j].parentTask(threadData[j].parentTaskData);
      threadData[j].parentIsFinished = true;
      
      signalCondition(threadData[j].parentTaskComplete);
    } 
    if (getNumElementsInQueue(&manager->threadQueue) != 0) {
      size_t j = pop(&manager->threadQueue);
      
      threadData[j].task = function;
      threadData[j].taskData = (data == NULL ? NULL : data[i]);
      manager->numThreadsRunning++;
      
      signalCondition(threadData[j].taskAvailable);
      
      ++i;
    }
  }
  
  while (manager->numThreadsRunning > 0) {
    waitOnCondition(manager->threadIsWaiting, manager->mutex);
    
    while (getNumElementsInQueue(&manager->parentTaskQueue) != 0) {
      size_t j = pop(&manager->parentTaskQueue);
      
      threadData[j].parentTask(threadData[j].parentTaskData);
      threadData[j].parentIsFinished = true;
      
      signalCondition(threadData[j].parentTaskComplete);
    }
  }
  
  unlockMutex(manager->mutex);
  
  return result;
}

int ext_btm_runTaskInParentThread(ext_btm_manager_t restrict manager, size_t threadId, ext_btm_taskFunction_t task, void* restrict data)
{
  if (manager->threads == NULL || manager->threadData == NULL ||
      manager->numThreadsActive == 0) return EINVAL;
  
  lockMutex(manager->mutex);
  
  ThreadData* threadData = manager->threadData + threadId;
  
  threadData->parentTask = task;
  threadData->parentTaskData = data;
  threadData->parentIsFinished = false;
      
  push(&manager->parentTaskQueue, threadId);
  
  signalCondition(manager->threadIsWaiting);
  
  while (!threadData->parentIsFinished) waitOnCondition(threadData->parentTaskComplete, manager->mutex);
  
  unlockMutex(manager->mutex);
  
  return 0;
}

static void* threadLoop(void* v_data)
{
  ThreadData* data = (ThreadData*) v_data;
  ext_btm_manager_t manager = data->manager;
  
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
    signalCondition(manager->threadIsWaiting);
    push(&manager->threadQueue, data->id);
  }
  
  manager->numThreadsActive--;
  
  unlockMutex(manager->mutex);
  
  return NULL;
}

int ext_btm_destroy(ext_btm_manager_t manager)
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
  
  result |= destroyCondition(manager->threadIsWaiting);
  result |= destroyCondition(manager->threadIsActive);
  result |= destroyMutex(manager->mutex);
  
  
  invalidateIndexArrayQueue(&manager->parentTaskQueue);
  invalidateIndexArrayQueue(&manager->threadQueue);
  
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  
  if (manager->threadData != NULL) {
    for (size_t i = 0; i < manager->numThreads; ++i) result |= destroyThreadData(&manager->threadData[i]);
    
    free(manager->threadData);
    manager->threadData = NULL;
  }
  
  free(manager);
  
  return result;
}

bool ext_btm_isNull(ext_btm_manager_t manager) {
  return manager == NULL;
}

static int initializeManager(ext_btm_manager_t manager, size_t numThreads)
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
  bool threadIsWaitingInitialized = false;
  
  manager->threads = (Thread*) malloc(numThreads * sizeof(Thread));
  if (manager->threads == NULL) { result = ENOMEM; goto ext_btm_manager_initialization_failed; }
  
  manager->threadData = (ThreadData*) malloc(numThreads * sizeof(ThreadData));
  if (manager->threadData == NULL) { result = ENOMEM; goto ext_btm_manager_initialization_failed; }
    
  result = initializeIndexArrayQueue(&manager->threadQueue, numThreads);
  if (result != 0) goto ext_btm_manager_initialization_failed;
  
  result = initializeIndexArrayQueue(&manager->parentTaskQueue, numThreads);
  if (result != 0) goto ext_btm_manager_initialization_failed;
  
  result = initializeMutex(manager->mutex);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) mutexInitialized = true;
    goto ext_btm_manager_initialization_failed;
  }
  
  result = initializeCondition(manager->threadIsActive);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) threadIsActiveInitialized = true;
    goto ext_btm_manager_initialization_failed;
  }
  
  result = initializeCondition(manager->threadIsWaiting);
  if (result != 0) {
    threadIsWaitingInitialized = true;
    goto ext_btm_manager_initialization_failed;
  }
  
  return result;
  
ext_btm_manager_initialization_failed:
  if (threadIsWaitingInitialized) destroyCondition(manager->threadIsWaiting);
  if (threadIsActiveInitialized) destroyCondition(manager->threadIsActive);
  if (mutexInitialized) destroyMutex(manager->mutex);
  
  invalidateIndexArrayQueue(&manager->parentTaskQueue);
  invalidateIndexArrayQueue(&manager->threadQueue);
  
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  if (manager->threadData != NULL) { free(manager->threadData); manager->threadData = NULL; }
  
  free(manager);
  
  return result;
}

static int initializeThreadData(ext_btm_manager_t manager, ThreadData* data, size_t id)
{
  data->manager = manager;
  data->id = id;
  
  data->task = NULL;
  data->taskData = NULL;
  
  data->parentTask = NULL;
  data->parentTaskData = NULL;
  data->parentIsFinished = true;
  
  bool taskAvailableInitialized = false;
  bool parentTaskCompleteInitialized = false;
  
  int result = initializeCondition(data->taskAvailable);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) taskAvailableInitialized = true;
    goto ext_btm_thread_initialization_failed;
  }
  
  result = initializeCondition(data->parentTaskComplete);
  if (result != 0) {
    parentTaskCompleteInitialized = true;
    goto ext_btm_thread_initialization_failed;
  }
  
  return result;
  
ext_btm_thread_initialization_failed:
  if (parentTaskCompleteInitialized) destroyCondition(data->parentTaskComplete);
  if (taskAvailableInitialized) destroyCondition(data->taskAvailable);
  
  return result;
}

static int destroyThreadData(ThreadData* data)
{
  return destroyCondition(data->taskAvailable);
}

/* static IndexArrayQueue* createIndexArrayQueue(size_t queueSize)
{
  IndexArrayQueue* result = (IndexArrayQueue*) malloc(sizeof(IndexArrayQueue));
  if (result == NULL) return result;
  
  if (initializeIndexArrayQueue(result, queueSize) != 0) {
    free(result);
    return NULL;
  }
  
  return result;
} */

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

/* static void destroyIndexArrayQueue(IndexArrayQueue* queue)
{
  invalidateIndexArrayQueue(queue);
  if (queue != NULL) free(queue);
} */


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

