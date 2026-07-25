#include "config.h"

#include <misc/thread.h>
#include "pthread.h"

#include "indexArrayQueue.h"     // IndexArrayQueue and its push/pop/etc.
#include "threadManagerCommon.h" // misc_partitionThreadJob

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>

struct ThreadData;

static int initializeManager(misc_btm_manager_t manager, size_t numThreads);
static int initializeThreadData(misc_btm_manager_t manager, struct ThreadData* data, size_t threadId);
static int destroyThreadData(struct ThreadData* data);

static void* threadLoop(void* _data);

typedef struct ThreadData {
  misc_btm_manager_t manager;
  Condition taskAvailable;
  Condition parentTaskComplete;
  size_t threadId;
  
  misc_btm_taskFunction_t task;
  void* taskData;
  
  misc_btm_taskFunction_t parentTask;
  void* parentTaskData;
  bool parentIsFinished;
} ThreadData;

typedef struct _misc_btm_manager_t {
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
  
} _misc_btm_manager_t;

int misc_btm_create(misc_btm_manager_t* managerPtr, size_t numThreads)
{
  *managerPtr = (misc_btm_manager_t) malloc(sizeof(_misc_btm_manager_t));
  if (*managerPtr == NULL) return ENOMEM;
  
  misc_btm_manager_t manager = *managerPtr;
  
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
    misc_btm_destroy(manager);
    *managerPtr = NULL;
  }
  
  return result;
}

size_t misc_btm_getThreadId(const misc_btm_manager_t manager)
{
  Thread nativeThreadId = pthread_self();
  size_t i;
  for (i = 0; i < manager->numThreads; ++i) if (nativeThreadId == manager->threads[i]) break;
  
  return i;
}

size_t misc_btm_getNumThreads(const misc_btm_manager_t manager)
{
  if (manager == NULL) return 1;
  return manager->numThreads;
}


void misc_btm_getNumThreadsForJob(const misc_btm_manager_t restrict threadManager, size_t numElements, size_t minNumElementsPerThread,
                                size_t* restrict numThreadsPtr, size_t* restrict numElementsPerThreadPtr, size_t* restrict offByOneIndexPtr)
{
  size_t numThreadsManaged = 0;
  if (numElements < 2 * minNumElementsPerThread || threadManager == NULL ||
      (numThreadsManaged = misc_btm_getNumThreads(threadManager)) <= 1) {
    if (numThreadsPtr != NULL) *numThreadsPtr = 1;
    *numElementsPerThreadPtr = numElements;
    *offByOneIndexPtr = 1;
    return;
  }

  misc_partitionThreadJob(numElements, minNumElementsPerThread, numThreadsManaged,
                          numThreadsPtr, numElementsPerThreadPtr, offByOneIndexPtr);
}

int misc_btm_runTasks(misc_btm_manager_t restrict manager, misc_btm_taskFunction_t function,
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

int misc_btm_runTaskInParentThread(misc_btm_manager_t restrict manager, size_t threadId, misc_btm_taskFunction_t task, void* restrict data)
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
  misc_btm_manager_t manager = data->manager;
  
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
    push(&manager->threadQueue, data->threadId);
  }
  
  manager->numThreadsActive--;
  
  unlockMutex(manager->mutex);
  
  return NULL;
}

int misc_btm_destroy(misc_btm_manager_t manager)
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

bool misc_btm_isNull(misc_btm_manager_t manager) {
  return manager == NULL;
}

static int initializeManager(misc_btm_manager_t manager, size_t numThreads)
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
  bool threadQueueInitialized = false;
  bool parentTaskQueueInitialized = false;
  
  manager->threads = (Thread*) malloc(numThreads * sizeof(Thread));
  if (manager->threads == NULL) { result = ENOMEM; goto misc_btm_manager_initialization_failed; }
  
  manager->threadData = (ThreadData*) malloc(numThreads * sizeof(ThreadData));
  if (manager->threadData == NULL) { result = ENOMEM; goto misc_btm_manager_initialization_failed; }
    
  result = initializeIndexArrayQueue(&manager->threadQueue, numThreads);
  if (result != 0) goto misc_btm_manager_initialization_failed;
  threadQueueInitialized = true;
  
  result = initializeIndexArrayQueue(&manager->parentTaskQueue, numThreads);
  if (result != 0) goto misc_btm_manager_initialization_failed;
  parentTaskQueueInitialized = true;
  
  result = initializeMutex(manager->mutex);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) mutexInitialized = true;
    goto misc_btm_manager_initialization_failed;
  }
  
  result = initializeCondition(manager->threadIsActive);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) threadIsActiveInitialized = true;
    goto misc_btm_manager_initialization_failed;
  }
  
  result = initializeCondition(manager->threadIsWaiting);
  if (result != 0) {
    threadIsWaitingInitialized = true;
    goto misc_btm_manager_initialization_failed;
  }
  
  return result;
  
misc_btm_manager_initialization_failed:
  if (threadIsWaitingInitialized) destroyCondition(manager->threadIsWaiting);
  if (threadIsActiveInitialized) destroyCondition(manager->threadIsActive);
  if (mutexInitialized) destroyMutex(manager->mutex);
  
  if (parentTaskQueueInitialized) invalidateIndexArrayQueue(&manager->parentTaskQueue);
  if (threadQueueInitialized) invalidateIndexArrayQueue(&manager->threadQueue);
  
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  if (manager->threadData != NULL) { free(manager->threadData); manager->threadData = NULL; }
  
  free(manager);
  
  return result;
}

static int initializeThreadData(misc_btm_manager_t manager, ThreadData* data, size_t threadId)
{
  data->manager = manager;
  data->threadId = threadId;
  
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
    goto misc_btm_thread_initialization_failed;
  }
  
  result = initializeCondition(data->parentTaskComplete);
  if (result != 0) {
    parentTaskCompleteInitialized = true;
    goto misc_btm_thread_initialization_failed;
  }
  
  return result;
  
misc_btm_thread_initialization_failed:
  if (parentTaskCompleteInitialized) destroyCondition(data->parentTaskComplete);
  if (taskAvailableInitialized) destroyCondition(data->taskAvailable);
  
  return result;
}

static int destroyThreadData(ThreadData* data)
{
  return destroyCondition(data->taskAvailable);
}

