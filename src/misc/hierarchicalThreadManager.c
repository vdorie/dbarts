#include "config.h"

#include <misc/thread.h>
#include "pthread.h"

#include <errno.h>
#include <stdarg.h> // varargs for buffered printf
#include <stdlib.h> // malloc

#include <external/io.h>

// clock_gettime + CLOCK_REALTIME are in time.h, gettimeofday is in sys/time.h; plain time() is in time.h too
// time.h imported from <misc/thread.h>
#if (!defined(HAVE_CLOCK_GETTIME) || !defined(CLOCK_REALTIME)) && defined(HAVE_GETTIMEOFDAY)
#  include <sys/time.h>
#endif

#ifdef __GNUC__
#  define UNUSED __attribute__ ((unused))
#else
#  define UNUSED
#endif

#define TASK_COMPLETE ((size_t) -1)
#define TASK_BEFORE_START (((size_t) -1) - 1)

#define BUFFER_LENGTH 32768

struct ThreadData;

typedef struct ThreadData {
  misc_htm_manager_t manager;
  size_t threadId;
  
  struct ThreadData* next;
  
  size_t topLevelTaskId;
  bool isTopLevelTask;
  
  union taskFunction_t {
    misc_htm_topLevelTaskFunction_t tl;
    misc_htm_subTaskFunction_t      st;
  } task;
  void* taskData;
    
  Condition taskAvailable;
} ThreadData;

typedef struct ThreadStack {
  ThreadData* first;
} ThreadStack;

typedef struct TopLevelTaskStatus {
  ThreadData* thread;
  ThreadStack threadStack;
  size_t numThreads; // num in stack + one base
  
  size_t progress;
  
  size_t numSubTaskPiecesInProgress;
  
  Condition taskDone;
} TopLevelTaskStatus;


typedef struct _misc_htm_manager_t {
  Thread* threads;
  size_t numThreads;
  
  ThreadData* threadData;
  
  TopLevelTaskStatus* topLevelTaskStatus;
  
  ThreadStack availableThreadStack;
  size_t numThreadsAvailable;
  
  size_t numTopLevelTasks;
  size_t numTopLevelTasksInProgress;
  
  Mutex mutex;
  Condition taskDone;
  
  char* buffer;
  size_t bufferPos;
  
  bool threadsShouldExit;
  
  Condition threadIsActive; // used to synchronize at start
} _misc_htm_manager_t;

static int initializeThreadData(misc_htm_manager_t manager, ThreadData* data, size_t threadId);
static int invalidateThreadData(ThreadData* thread);

static int initializeTopLevelTaskStatus(TopLevelTaskStatus* status);
static int invalidateTopLevelTaskStatus(TopLevelTaskStatus* status);

static void initializeThreadStack(ThreadStack* stack);
static void invalidateThreadStack(ThreadStack* stack);
static ThreadData* pop(ThreadStack* stack);
static ThreadData* popN(ThreadStack* stack, size_t n);
static void push(ThreadStack* stack, ThreadData* thread);
static size_t pushN(ThreadStack* stack, ThreadData* thread);
static bool stackIsEmpty(ThreadStack* stack);

UNUSED static void printManagerStatus(const misc_htm_manager_t manager);
static void printTaskProgress(TopLevelTaskStatus* status);


int misc_htm_runTopLevelTasks(misc_htm_manager_t restrict manager, misc_htm_topLevelTaskFunction_t function,
                             void** restrict data, size_t numTasks)
{
  if (manager->threads == NULL || manager->threadData == NULL) return EINVAL;
  
  lockMutex(manager->mutex);
  
  manager->topLevelTaskStatus = (TopLevelTaskStatus*) malloc(numTasks * sizeof(TopLevelTaskStatus));
  if (manager->topLevelTaskStatus == NULL) { unlockMutex(manager->mutex); return ENOMEM; }
  manager->numTopLevelTasks = numTasks;
  
  size_t taskId = 0;
  int result = 0;
  
  for (taskId = 0; taskId < numTasks; ++taskId)
    if ((result = initializeTopLevelTaskStatus(&manager->topLevelTaskStatus[taskId])) != 0) break;
  
  if (result != 0) {
    for ( /* */ ; taskId > 0; --taskId)
      invalidateTopLevelTaskStatus(&manager->topLevelTaskStatus[taskId - 1]);
    free(manager->topLevelTaskStatus);
    unlockMutex(manager->mutex);
    
    return result;
  }
  
  // ext_printf("running top %lu level tasks without output\n", numTasks);
  
  for (taskId = 0; taskId < numTasks; ++taskId) {
    while (stackIsEmpty(&manager->availableThreadStack))
      waitOnCondition(manager->taskDone, manager->mutex);
    
    ThreadData* thread = pop(&manager->availableThreadStack);
    manager->numThreadsAvailable--;
    
    manager->topLevelTaskStatus[taskId].thread = thread;
    manager->topLevelTaskStatus[taskId].numThreads = 1;
    
    thread->task.tl = function;
    thread->taskData = (data == NULL ? NULL : data[taskId]);
    thread->topLevelTaskId = taskId;
    thread->isTopLevelTask = true;
    
    manager->numTopLevelTasksInProgress++;
    
    signalCondition(thread->taskAvailable);
  }
  
  // ext_printf("waiting for top level tasks to finish\n");
  // printManagerStatus(manager);
  
  while (manager->numTopLevelTasksInProgress > 0)
    waitOnCondition(manager->taskDone, manager->mutex);
  
  // ext_printf("cleaning up top level tasks\n");
  
  for ( /* */ ; taskId > 0; --taskId)
    result |= invalidateTopLevelTaskStatus(&manager->topLevelTaskStatus[taskId - 1]);
  free(manager->topLevelTaskStatus);
  manager->topLevelTaskStatus = NULL;
  manager->numTopLevelTasks = 0;
  
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

int misc_htm_runTopLevelTasksWithOutput(misc_htm_manager_t restrict manager, misc_htm_topLevelTaskFunction_t function,
                                       void** restrict data, size_t numTasks, const struct timespec* restrict outputDelay)
{
  if (manager->threads == NULL || manager->threadData == NULL) return EINVAL;
  
  struct timespec wakeTime;
  
  lockMutex(manager->mutex);
  
  manager->topLevelTaskStatus = (TopLevelTaskStatus*) malloc(numTasks * sizeof(TopLevelTaskStatus));
  if (manager->topLevelTaskStatus == NULL) { unlockMutex(manager->mutex); return ENOMEM; }
  manager->numTopLevelTasks = numTasks;
  
  size_t taskId = 0;
  int result = 0;
  
  for (taskId = 0; taskId < numTasks; ++taskId)
    if ((result = initializeTopLevelTaskStatus(&manager->topLevelTaskStatus[taskId])) != 0) break;
  
  if (result != 0) {
    for ( /* */ ; taskId > 0; --taskId)
      invalidateTopLevelTaskStatus(&manager->topLevelTaskStatus[taskId - 1]);
    free(manager->topLevelTaskStatus);
    unlockMutex(manager->mutex);
    
    return result;
  }
  
  getTime(&wakeTime);
  // ext_printf("wake time: %ld.%.9ld", (long long int) wakeTime.tv_sec, (long long int) wakeTime.tv_nsec);
  wakeTime.tv_sec  += outputDelay->tv_sec;
  wakeTime.tv_nsec += outputDelay->tv_nsec;
  // ext_printf(", added time: %ld.%.9ld", (long long int) wakeTime.tv_sec, (long long int) wakeTime.tv_nsec);
  wakeTime.tv_sec  += wakeTime.tv_nsec / 1000000000;
  wakeTime.tv_nsec %= 1000000000;
  // ext_printf(", mod time: %ld.%.9ld\n", (long long int) wakeTime.tv_sec, (long long int) wakeTime.tv_nsec);
  
  // ext_printf("running top %lu level tasks with output\n", numTasks);
  
  for (taskId = 0; taskId < numTasks; ++taskId) {
    // while (stackIsEmpty(&manager->availableThreadStack)) waitOnCondition(manager->taskDone, manager->mutex);
    while (stackIsEmpty(&manager->availableThreadStack)) {
      int waitStatus = waitOnConditionForTime(manager->taskDone, manager->mutex, wakeTime);
      if (waitStatus == ETIMEDOUT) {
        if (manager->bufferPos != 0) {
          ext_printf(manager->buffer);
          ext_fflush_stdout();
          manager->bufferPos = 0;
        }
        
        getTime(&wakeTime);
        wakeTime.tv_sec  += outputDelay->tv_sec;
        wakeTime.tv_nsec += outputDelay->tv_nsec;
        wakeTime.tv_sec  += wakeTime.tv_nsec / 1000000000;
        wakeTime.tv_nsec %= 1000000000;
      }
    }
    
    ThreadData* thread = pop(&manager->availableThreadStack);
    manager->numThreadsAvailable--;
    
    manager->topLevelTaskStatus[taskId].thread = thread;
    manager->topLevelTaskStatus[taskId].numThreads = 1;
    
    thread->task.tl = function;
    thread->taskData = (data == NULL ? NULL : data[taskId]);
    thread->topLevelTaskId = taskId;
    thread->isTopLevelTask = true;
    
    manager->numTopLevelTasksInProgress++;
    
    signalCondition(thread->taskAvailable);
  }
  
  // ext_printf("waiting for top level tasks to finish\n");
  // printManagerStatus(manager);
  
  while (manager->numTopLevelTasksInProgress > 0) {
    int waitStatus = waitOnConditionForTime(manager->taskDone, manager->mutex, wakeTime);
    if (waitStatus == ETIMEDOUT) {
      if (manager->bufferPos != 0) {
        ext_printf(manager->buffer);
        ext_fflush_stdout();
        manager->bufferPos = 0;
      }
      
      getTime(&wakeTime);
      wakeTime.tv_sec  += outputDelay->tv_sec;
      wakeTime.tv_nsec += outputDelay->tv_nsec;
      wakeTime.tv_sec  += wakeTime.tv_nsec / 1000000000;
      wakeTime.tv_nsec %= 1000000000;
    }
  }
  
  // ext_printf("cleaning up top level tasks\n");
  
  for ( /* */ ; taskId > 0; --taskId)
    result |= invalidateTopLevelTaskStatus(&manager->topLevelTaskStatus[taskId - 1]);
  free(manager->topLevelTaskStatus);
  manager->topLevelTaskStatus = NULL;
  manager->numTopLevelTasks = 0;
  
  if (manager->bufferPos != 0) {
    ext_printf(manager->buffer);
    ext_fflush_stdout();
    manager->bufferPos = 0;
  }
  
  unlockMutex(manager->mutex);
  
  return result;
}


int misc_htm_runSubTask(misc_htm_manager_t restrict manager, size_t taskId, misc_htm_subTaskFunction_t subTask,
                       void** restrict data, size_t numPieces)
{
  if (manager->threads == NULL || manager->threadData == NULL || manager->topLevelTaskStatus == NULL) return EINVAL;
  
  TopLevelTaskStatus* taskStatus = &manager->topLevelTaskStatus[taskId];
  ThreadData* thread = taskStatus->threadStack.first;
  
  lockMutex(manager->mutex);
  
  if (numPieces > 1) for (size_t i = 1; i < numPieces; ++i) {
    
    thread->task.st = subTask;
    thread->taskData = (data == NULL ? NULL : data[i]);
    thread->topLevelTaskId = taskId;
    thread->isTopLevelTask = false;
    
    taskStatus->numSubTaskPiecesInProgress++;
    
    signalCondition(thread->taskAvailable);
    
    thread = thread->next;
  }
  
  unlockMutex(manager->mutex);
  
  subTask(data[0]);
  
  lockMutex(manager->mutex);
  
  while (taskStatus->numSubTaskPiecesInProgress > 0) waitOnCondition(taskStatus->taskDone, manager->mutex);
  
  unlockMutex(manager->mutex);
  
  return 0;
}

misc_size_t misc_htm_getNumThreadsForTopLevelTask(const misc_htm_manager_t threadManager, size_t taskId)
{
  if (threadManager == NULL || taskId == EXT_HTM_INVALID_TASK_ID || threadManager->topLevelTaskStatus == NULL) return 1;
  return threadManager->topLevelTaskStatus[taskId].numThreads;
}

void misc_htm_getNumPiecesForSubTask(const misc_htm_manager_t restrict threadManager, size_t taskId,
                                    size_t numElements, size_t minNumElementsPerPiece,
                                    size_t* restrict numPiecesPtr, size_t* restrict numElementsPerPiecePtr, size_t* restrict offByOneIndexPtr)
{
  size_t numThreadsManaged = 0;
  if (numElements < 2 * minNumElementsPerPiece || threadManager == NULL || taskId == EXT_HTM_INVALID_TASK_ID ||
      (numThreadsManaged = threadManager->topLevelTaskStatus[taskId].numThreads) <= 1) {
    if (numPiecesPtr != NULL) *numPiecesPtr = 1;
    *numElementsPerPiecePtr = numElements;
    *offByOneIndexPtr = 1;
    return;
  }
  
  size_t numPieces = minNumElementsPerPiece == 0 ? numElements : numElements / minNumElementsPerPiece;
  if (numPieces > numThreadsManaged) numPieces = numThreadsManaged;
  
  size_t numElementsPerPiece = numElements / numPieces;
  size_t offByOneIndex = numElements % numPieces;
  if (offByOneIndex != 0) numElementsPerPiece += 1;
  else offByOneIndex = numPieces;
  
  if (numPiecesPtr != NULL) *numPiecesPtr = numPieces;
  *numElementsPerPiecePtr = numElementsPerPiece;
  *offByOneIndexPtr = offByOneIndex;
}

static void* threadLoop(void* _thread)
{
  ThreadData* thread = (ThreadData*) _thread;
  misc_htm_manager_t manager = thread->manager;
  
  lockMutex(manager->mutex);
  
  push(&manager->availableThreadStack, thread);
  manager->numThreadsAvailable++;
  
  signalCondition(manager->threadIsActive);
  
  while (true) {
    while (thread->task.tl == NULL && manager->threadsShouldExit == false)
      waitOnCondition(thread->taskAvailable, manager->mutex);
    if (manager->threadsShouldExit == true) break;
    
    unlockMutex(manager->mutex);
    
    if (thread->isTopLevelTask)
      thread->task.tl(thread->topLevelTaskId, thread->taskData);
    else
      thread->task.st(thread->taskData);
    
    lockMutex(manager->mutex);
    
    thread->task.tl = NULL;
    thread->taskData = NULL;
    TopLevelTaskStatus* taskStatus = &manager->topLevelTaskStatus[thread->topLevelTaskId];
    
    if (thread->isTopLevelTask) {
      // return both the main thread and all subthreads to the manager
      thread->next = popN(&taskStatus->threadStack, taskStatus->numThreads - 1);
      
      pushN(&manager->availableThreadStack, thread);
      
      manager->numThreadsAvailable += taskStatus->numThreads;
      
      manager->numTopLevelTasksInProgress--;
      taskStatus->progress = TASK_COMPLETE;
      taskStatus->thread = NULL;
      
      // ext_printf("task %lu complete, popping %lu threads\n", thread->topLevelTaskId, taskStatus->numThreads - 1);
      // printManagerStatus(manager);
      
      signalCondition(manager->taskDone);
    } else {
      taskStatus->numSubTaskPiecesInProgress--;
      signalCondition(taskStatus->taskDone);
    }
  }
  
  unlockMutex(manager->mutex);
  
  return NULL;
}

misc_size_t misc_htm_reserveThreadsForSubTask(misc_htm_manager_t manager, size_t taskId, size_t progress)
{
  size_t numTasksMoreComplete = 0;
  
  lockMutex(manager->mutex);
  
  // count how many threads this task should get based on all tasks' progress
  for (size_t i = 0; i < manager->numTopLevelTasks; ++i) {
    if (manager->topLevelTaskStatus[i].progress >= progress &&
        manager->topLevelTaskStatus[i].progress != TASK_COMPLETE &&
        manager->topLevelTaskStatus[i].progress != TASK_BEFORE_START)
    {
      ++numTasksMoreComplete;
    }
  }
  manager->topLevelTaskStatus[taskId].progress = progress;
  
  size_t minorShare = manager->numThreads / manager->numTopLevelTasksInProgress;
  size_t majorShare = minorShare + 1;
  size_t numWithMinorShare = majorShare * manager->numTopLevelTasks - manager->numThreads;
  
  size_t newNumThreads = numTasksMoreComplete < numWithMinorShare ? minorShare : majorShare;
  
  TopLevelTaskStatus* taskStatus = &manager->topLevelTaskStatus[taskId];
  
  if (newNumThreads == taskStatus->numThreads) {
    // ext_printf("task id %lu: no change\n", taskId);
    unlockMutex(manager->mutex);
    return newNumThreads;
  }
  
  // ext_printf("task id %lu: %lu -> %lu threads\n", taskId, taskStatus->numThreads, newNumThreads);
  
  if (newNumThreads > taskStatus->numThreads) {
    
    ThreadData* threads = popN(&manager->availableThreadStack, newNumThreads - taskStatus->numThreads);
    
    pushN(&taskStatus->threadStack, threads);
    
    manager->numThreadsAvailable -= newNumThreads - taskStatus->numThreads;
    taskStatus->numThreads = newNumThreads;
    
  } else {
    
    ThreadData* threads = popN(&taskStatus->threadStack, taskStatus->numThreads - newNumThreads);
    
    pushN(&manager->availableThreadStack, threads);
        
    manager->numThreadsAvailable += taskStatus->numThreads - newNumThreads;
    taskStatus->numThreads = newNumThreads;
  }
  
  // printManagerStatus(manager);
  
  unlockMutex(manager->mutex);
  
  return newNumThreads;
}

static int initializeManager(misc_htm_manager_t manager, size_t numThreads)
{
  int result;
  
  manager->threads = NULL;
  manager->numThreads = numThreads;
  
  manager->threadData = NULL;
  manager->topLevelTaskStatus = NULL;
  
  initializeThreadStack(&manager->availableThreadStack);
  manager->numThreadsAvailable = 0;
  
  manager->numTopLevelTasks = 0;
  manager->numTopLevelTasksInProgress = 0;
  
  manager->threadsShouldExit = false;
  
  manager->buffer = NULL;
  
  bool mutexInitialized = false;
  bool threadIsActiveInitialized = false;
  bool taskDoneInitialized = false;
  
  manager->threads = (Thread*) malloc(numThreads * sizeof(Thread));
  if (manager->threads == NULL) { result = ENOMEM; goto misc_htm_initialization_failed; }
  
  manager->threadData = (ThreadData*) malloc(numThreads * sizeof(ThreadData));
  if (manager->threadData == NULL) { result = ENOMEM; goto misc_htm_initialization_failed; }
  
  manager->buffer = (char*) malloc(BUFFER_LENGTH * sizeof(char));
  if (manager->buffer == NULL) { result = ENOMEM; goto misc_htm_initialization_failed; }
  manager->bufferPos = 0;
      
  result = initializeMutex(manager->mutex);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) mutexInitialized = true;
    goto misc_htm_initialization_failed;
  }
  
  result = initializeCondition(manager->threadIsActive);
  if (result != 0) {
    if (result != EBUSY && result != EINVAL) threadIsActiveInitialized = true;
    goto misc_htm_initialization_failed;
  }
  
  result = initializeCondition(manager->taskDone);
  if (result != 0) {
    taskDoneInitialized = true;
    goto misc_htm_initialization_failed;
  }
  
  return result;
  
misc_htm_initialization_failed:
  if (manager->buffer != NULL) { free(manager->buffer); manager->buffer = NULL; }
  if (manager->threadData != NULL) { free(manager->threadData); manager->threadData = NULL; }
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  
  // invalidateThreadStack(&manager->availableThreadStack); // not needed at the moment
  
  if (mutexInitialized) destroyMutex(manager->mutex);
  if (threadIsActiveInitialized) destroyCondition(manager->threadIsActive);
  if (taskDoneInitialized) destroyCondition(manager->taskDone);
  
  free(manager);
  
  return result;
}

int misc_htm_create(misc_htm_manager_t* managerPtr, size_t numThreads)
{
  *managerPtr = (misc_htm_manager_t) malloc(sizeof(_misc_htm_manager_t));
  if (*managerPtr == NULL) return ENOMEM;
  
  misc_htm_manager_t manager = *managerPtr;
  
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
      invalidateThreadData(&manager->threadData[i]);
      manager->numThreads = i;
      break;
    }
  }
  
  // wait for all threads to check in
  lockMutex(manager->mutex);
  while (manager->numThreadsAvailable < manager->numThreads)
    waitOnCondition(manager->threadIsActive, manager->mutex);
  unlockMutex(manager->mutex);
  
  if (result != 0) {
    misc_htm_destroy(manager);
    *managerPtr = NULL;
  }
  
  return result;
}

int misc_htm_destroy(misc_htm_manager_t manager)
{
  if (manager == NULL) return 0;
    
  int result = 0;
  
  if (manager->topLevelTaskStatus != NULL && manager->numTopLevelTasks > 0) {
    // this really shouldn't be happening
    lockMutex(manager->mutex);
    while (manager->numTopLevelTasksInProgress > 0) waitOnCondition(manager->taskDone, manager->mutex);
  
    for ( size_t taskId = manager->numTopLevelTasks; taskId > 0; --taskId)
      result |= invalidateTopLevelTaskStatus(&manager->topLevelTaskStatus[taskId - 1]);
    free(manager->topLevelTaskStatus);
    manager->topLevelTaskStatus = NULL;
    manager->numTopLevelTasks = 0;
    
    unlockMutex(manager->mutex);
  }
  
  // check that initialization didn't fail
  if (manager->threads != NULL && manager->threadData != NULL &&
      manager->numThreadsAvailable > 0 && manager->numThreads > 0)
  {
    lockMutex(manager->mutex);
    
    manager->threadsShouldExit = true;
    
    for (size_t i = 0; i < manager->numThreads; ++i)
      signalCondition(manager->threadData[i].taskAvailable);
    
    unlockMutex(manager->mutex);
    
    for (size_t i = 0; i < manager->numThreads; ++i)
      result |= joinThread(manager->threads[i]);
  }
  
  result |= destroyCondition(manager->taskDone);
  result |= destroyCondition(manager->threadIsActive);
  result |= destroyMutex(manager->mutex);
  
  if (manager->buffer != NULL) {
    free(manager->buffer);
    manager->buffer = NULL;
  }
  
  
  if (manager->threadData != NULL) {
    for (size_t i = 0; i < manager->numThreads; ++i) result |= invalidateThreadData(&manager->threadData[i]);
    
    free(manager->threadData);
    manager->threadData = NULL;
  }
  
  if (manager->threads != NULL) { free(manager->threads); manager->threads = NULL; }
  
  invalidateThreadStack(&manager->availableThreadStack);
  
  free(manager);
  
  return result;
}

static int initializeThreadData(misc_htm_manager_t manager, ThreadData* data, size_t threadId)
{
  data->manager = manager;
  data->threadId = threadId;
  
  data->next = NULL;
  data->topLevelTaskId = EXT_HTM_INVALID_TASK_ID;
  data->isTopLevelTask = false;
  
  data->task.tl = NULL;
  data->taskData = NULL;
  
  int result = initializeCondition(data->taskAvailable);
  
  if (result != 0 && result != EBUSY && result != EINVAL)
    destroyCondition(data->taskAvailable);
  
  return result;
}

static int invalidateThreadData(ThreadData* data)
{
  return destroyCondition(data->taskAvailable);
}

static int initializeTopLevelTaskStatus(TopLevelTaskStatus* status)
{
  status->thread = NULL;
  
  status->numThreads = 0;
  status->progress = TASK_BEFORE_START;
  status->numSubTaskPiecesInProgress = 0;
  status->threadStack.first = NULL;
  
  int result = initializeCondition(status->taskDone);
  
  if (result != 0 && result != EBUSY && result != EINVAL)
    destroyCondition(status->taskDone);
  
  return result;
}

static int invalidateTopLevelTaskStatus(TopLevelTaskStatus* status)
{
  return destroyCondition(status->taskDone);
}

static void initializeThreadStack(ThreadStack* stack)
{
  stack->first = NULL;
}

static void invalidateThreadStack(ThreadStack* stack)
{
  stack->first = NULL;
}

static ThreadData* pop(ThreadStack* stack)
{
  ThreadData* result = stack->first;
  
  stack->first = result->next;
  result->next = NULL;
  
  return result;
}

static ThreadData* popN(ThreadStack* stack, size_t n)
{
  if (n == 0) return NULL;
  
  ThreadData* first = stack->first;
  ThreadData* last = first;
  
  for (size_t i = 1; i < n; ++i)
    last = last->next;
  
  stack->first = last->next;
  last->next = NULL;
  
  return first;
}

static void push(ThreadStack* stack, ThreadData* thread)
{
  thread->next = stack->first;
  stack->first = thread;
}

static size_t pushN(ThreadStack* stack, ThreadData* thread)
{
  if (thread == NULL) return 0;
  
  ThreadData* last = thread;
  size_t n = 1;
  
  while (last->next != NULL) {
    last = last->next;
    ++n;
  }
  
  last->next = stack->first;
  stack->first = thread;
  
  return n;
}

static bool stackIsEmpty(ThreadStack* stack) {
  return stack->first == NULL;
}

void misc_htm_printf(misc_htm_manager_t manager, const char* format, ...)
{
  if (manager == NULL || manager->numThreads == 0) {
    char buffer[BUFFER_LENGTH];
    
    va_list argsPointer;
    va_start(argsPointer, format);
    vsnprintf(buffer, BUFFER_LENGTH, format, argsPointer);
    va_end(argsPointer);
    
    ext_printf(buffer);
    
    return;
  }
  
  lockMutex(manager->mutex);
  
  va_list argsPointer;
  va_start(argsPointer, format);
  manager->bufferPos += vsnprintf(manager->buffer + manager->bufferPos, BUFFER_LENGTH - manager->bufferPos, format, argsPointer);
  va_end(argsPointer);
  
  unlockMutex(manager->mutex);
}

// debug functions...
static void printThreadStack(ThreadStack* stack) {
  if (stack->first == NULL) {
    ext_printf("empty");
    return;
  }
  
  ext_printf("%lu", stack->first->threadId);
  
  ThreadData* thread = stack->first->next;
  while (thread != NULL) {
    ext_printf(", %lu", thread->threadId);
    thread = thread->next;
  }
}

static void printManagerStatus(const misc_htm_manager_t manager)
{
  if (manager->numTopLevelTasks == 0) {
    ext_printf("status: inactive\n");
    ext_printf("  avail pool: ");
    printThreadStack(&manager->availableThreadStack);
    ext_printf("\n\n");
  } else {
    ext_printf("status: running %lu tasks, %lu in progress\n", manager->numTopLevelTasks, manager->numTopLevelTasksInProgress);
    for (size_t i = 0; i < manager->numTopLevelTasks; ++i) {
      TopLevelTaskStatus* taskStatus = &manager->topLevelTaskStatus[i];
      
      ext_printf("  task %lu: ", i);
      printTaskProgress(taskStatus);
      ext_printf(" progress, %lu subtasks, threads ", taskStatus->numSubTaskPiecesInProgress);
      
      if (taskStatus->thread == NULL) ext_printf("none");
      else ext_printf("%lu", taskStatus->thread->threadId);
      
      ThreadData* thread = taskStatus->threadStack.first;
      while (thread != NULL) {
        ext_printf(" %lu", thread->threadId);
        thread = thread->next;
      }
      ext_printf("\n");
    }
    
    ext_printf("  avail pool: ");
    printThreadStack(&manager->availableThreadStack);
    ext_printf("\n\n");
  }
}

static void printTaskProgress(TopLevelTaskStatus* status)
{
  if (status->progress == TASK_COMPLETE) ext_printf("complete");
  else if (status->progress == TASK_BEFORE_START) ext_printf("before start");
  else ext_printf("%lu", status->progress);
}

