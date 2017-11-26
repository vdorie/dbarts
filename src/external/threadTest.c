#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <external/thread.h>

#define NUM_THREADS 4
#define NUM_TASKS   3

#define NUM_REPETITIONS 6
// #define NUM_ITERATIONS 100000000

#define FIND_PRIME 100000lu

static uint64_t findPrimeNumber(uint64_t n)
{
  uint64_t count = 0;
  uint64_t a = 2;
  
  while (count < n) {
    uint64_t b = 2;
    bool prime = true;
    
    while (b * b <= a) {
      if (a % b == 0) {
        prime = false;
        break;
      }
      ++b;
    }
    if (prime) ++count;
    ++a;
  }
  return --a;
}

typedef struct TaskData {
  ext_htm_manager_t manager;
} TaskData;

typedef struct SubTaskData {
  size_t numThreads;
} SubTaskData;

static void subTaskFunction(void* _data) {
  SubTaskData* subTaskData = (SubTaskData*) _data;
  
  size_t numPrimes = FIND_PRIME / subTaskData->numThreads;
  
  findPrimeNumber(numPrimes);
}

static void topLevelTaskFunction(size_t taskId, void* _data) {
  TaskData* data = (TaskData*) _data;
  
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    size_t numThreads = ext_htm_reserveThreadsForSubTask(data->manager, taskId, i);
    
    SubTaskData* subTaskData = (SubTaskData*) malloc(numThreads * sizeof(SubTaskData));
    void** subTaskDataPtr = (void**) malloc(numThreads * sizeof(void*));
    
    for (size_t j = 0; j < numThreads; ++j) {
      subTaskData[j].numThreads = numThreads;
      subTaskDataPtr[j] = &subTaskData[j];
    }
    
    ext_htm_runSubTasks(data->manager, &subTaskFunction, subTaskDataPtr, taskId);
    
    free(subTaskDataPtr);
    free(subTaskData);
  }
}



int main(int argc, const char* argv[]) {
  ext_htm_manager_t manager;
  
  int result = ext_htm_create(&manager, NUM_THREADS);
  if (result != 0) {
    printf("error %s", strerror(result));
    return result;
  }
  
  TaskData* taskData = (TaskData*) malloc(NUM_TASKS * sizeof(TaskData));
  void** taskDataPtrs = (void**) malloc(NUM_TASKS * sizeof(void*));
  for (size_t i = 0; i < NUM_TASKS; ++i) {
    taskData[i].manager = manager;
    taskDataPtrs[i] = &taskData[i];
  }
  
  
  ext_htm_runTopLevelTasks(manager, &topLevelTaskFunction, taskDataPtrs, NUM_TASKS);
  
  free(taskDataPtrs);
  free(taskData);
  
  ext_htm_destroy(manager);
}

