#include "config.h"
#include <misc/thread.h>
#include "pthread.h"

#include "indexArrayQueue.h"     // IndexArrayQueue and its push/pop/etc.
#include "threadManagerCommon.h" // misc_partitionThreadJob

#include <stdlib.h>
#include <errno.h>
#include <stdbool.h>

// clock_gettime + CLOCK_REALTIME are in time.h, gettimeofday is in sys/time.h; plain time() is in time.h too
// time.h imported from <misc/thread.h>
#if (!defined(HAVE_CLOCK_GETTIME) || !defined(CLOCK_REALTIME)) && defined(HAVE_GETTIMEOFDAY)
#  include <sys/time.h>
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

struct ThreadData;

static int initializeManager(misc_mt_manager_t manager, size_t numThreads);
static int initializeThreadData(misc_mt_manager_t manager, struct ThreadData* data, size_t threadId);
static int destroyThreadData(struct ThreadData* data);

static void* threadLoop(void* _data);

typedef struct ThreadData {
  misc_mt_manager_t manager;
  Condition taskAvailable;
  size_t threadId;
  
  misc_mt_taskFunction_t task;
  void* taskData;
} ThreadData;

typedef struct _misc_mt_manager_t {
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
} _misc_mt_manager_t;


int misc_mt_create(misc_mt_manager_t* managerPtr, size_t numThreads)
{
  *managerPtr = (misc_mt_manager_t) malloc(sizeof(_misc_mt_manager_t));
  if (*managerPtr == NULL) return ENOMEM;
  
  misc_mt_manager_t manager = *managerPtr;
  
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
    misc_mt_destroy(manager);
    *managerPtr = NULL;
  }
  
  return result;
}

size_t misc_mt_getNumThreads(const misc_mt_manager_t manager)
{
  if (manager == NULL) return 1;
  return manager->numThreads;
}


void misc_mt_getNumThreadsForJob(const misc_mt_manager_t restrict threadManager, size_t numElements, size_t minNumElementsPerThread,
                                 size_t* restrict numThreadsPtr, size_t* restrict numElementsPerThreadPtr, size_t* restrict offByOneIndexPtr)
{
  size_t numThreadsManaged = 0;
  if (numElements < 2 * minNumElementsPerThread || threadManager == NULL ||
      (numThreadsManaged = misc_mt_getNumThreads(threadManager)) <= 1) {
    if (numThreadsPtr != NULL) *numThreadsPtr = 1;
    *numElementsPerThreadPtr = numElements;
    *offByOneIndexPtr = 1;
    return;
  }

  misc_partitionThreadJob(numElements, minNumElementsPerThread, numThreadsManaged,
                          numThreadsPtr, numElementsPerThreadPtr, offByOneIndexPtr);
}

int misc_mt_runTasks(misc_mt_manager_t restrict manager, misc_mt_taskFunction_t function,
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

static inline void getTime(struct timespec* ts)
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

int misc_mt_runTasksWithInfo(misc_mt_manager_t restrict manager, misc_mt_taskFunction_t function,
                            void** restrict data, size_t numTasks, time_t sleepSeconds, misc_mt_infoFunction_t info)
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
  misc_mt_manager_t manager = data->manager;
  
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
    push(&manager->threadQueue, data->threadId);
  }
  
  manager->numThreadsActive--;
  
  unlockMutex(manager->mutex);
  
  return NULL;
}

int misc_mt_destroy(misc_mt_manager_t manager)
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

static int initializeManager(misc_mt_manager_t manager, size_t numThreads)
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
  bool threadQueueInitialized = false;
  
  manager->threads = (Thread*) malloc(numThreads * sizeof(Thread));
  if (manager->threads == NULL) { result = ENOMEM; goto misc_mt_initialization_failed; }
  
  manager->threadData = (ThreadData*) malloc(numThreads * sizeof(ThreadData));
  if (manager->threadData == NULL) { result = ENOMEM; goto misc_mt_initialization_failed; }
    
  result = initializeIndexArrayQueue(&manager->threadQueue, numThreads);
  if (result != 0) goto misc_mt_initialization_failed;
  threadQueueInitialized = true;
  
  result = initializeMutex(manager->mutex);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) mutexInitialized = true;
    goto misc_mt_initialization_failed;
  }
  
  result = initializeCondition(manager->threadIsActive);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) threadIsActiveInitialized = true;
    goto misc_mt_initialization_failed;
  }
  
  result = initializeCondition(manager->taskDone);
  if (result != 0) {
    taskDoneInitialized = true;
    goto misc_mt_initialization_failed;
  }
  
  return result;
  
misc_mt_initialization_failed:
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  if (manager->threadData != NULL) { free(manager->threadData); manager->threadData = NULL; }
  
  if (threadQueueInitialized) invalidateIndexArrayQueue(&manager->threadQueue);
  
  if (mutexInitialized) destroyMutex(manager->mutex);
  if (threadIsActiveInitialized) destroyCondition(manager->threadIsActive);
  if (taskDoneInitialized) destroyCondition(manager->taskDone);
  
  free(manager);
  
  return result;
}

static int initializeThreadData(misc_mt_manager_t manager, ThreadData* data, size_t threadId)
{
  data->manager = manager;
  data->threadId = threadId;
  
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
