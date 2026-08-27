#ifndef MISC_THREAD_H
#define MISC_THREAD_H

#ifdef __INTEL_COMPILER
#  define __need_timespec 1
#endif

#include <stdbool.h>
#include <misc/stddef.h>
#include <sys/types.h> // time_t
#include <time.h> // struct timespec

#ifdef __cplusplus
extern "C" {
#endif

struct _misc_mt_manager_t;
typedef struct _misc_mt_manager_t* misc_mt_manager_t;


// on failure, sets manager to NULL and returns an error code
int misc_mt_create(misc_mt_manager_t* manager, misc_size_t numThreads);
int misc_mt_destroy(misc_mt_manager_t manager);

misc_size_t misc_mt_getNumThreads(const misc_mt_manager_t manager);
// assign numElementsPerThread for the first offByOneIndex threads, and numElementsPerThread - 1 for those after
void misc_mt_getNumThreadsForJob(const misc_mt_manager_t restrict threadManager, misc_size_t numElements, misc_size_t minNumElementsPerThread,
                                 misc_size_t* restrict numThreadsPtr, misc_size_t* restrict numElementsPerThreadPtr, misc_size_t* restrict offByOneIndexPtr);


typedef void (*misc_mt_taskFunction_t)(void*);

typedef void (*misc_mt_infoFunction_t)(void** threadData, misc_size_t numThreads);

int misc_mt_runTasks(misc_mt_manager_t restrict manager, misc_mt_taskFunction_t task,
                     void** restrict data, misc_size_t numTasks);

int misc_mt_runTasksWithInfo(misc_mt_manager_t restrict manager, misc_mt_taskFunction_t task,
                             void** restrict data, misc_size_t numTasks, time_t sleepSeconds,
                             misc_mt_infoFunction_t info);


#ifdef __cplusplus
}
#endif

#endif // MISC_THREAD_H

