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


// hierarchical thread manager
#define EXT_HTM_INVALID_TASK_ID ((misc_size_t) -1)
struct _misc_htm_manager_t;
typedef struct _misc_htm_manager_t* misc_htm_manager_t;


typedef void (*misc_htm_topLevelTaskFunction_t)(misc_size_t taskId, void*);
typedef void (*misc_htm_subTaskFunction_t)(void*);

int misc_htm_create(misc_htm_manager_t* manager, misc_size_t numThreads);
int misc_htm_destroy(misc_htm_manager_t manager);

misc_size_t misc_htm_getNumThreads(const misc_htm_manager_t manager);

int misc_htm_runTopLevelTasks(misc_htm_manager_t restrict manager, misc_htm_topLevelTaskFunction_t task,
                              void** restrict data, misc_size_t numTasks);
int misc_htm_runTopLevelTasksWithOutput(misc_htm_manager_t restrict manager, misc_htm_topLevelTaskFunction_t task,
                                        void** restrict data, misc_size_t numTasks, const struct timespec* restrict outputDelay);

misc_size_t misc_htm_reserveThreadsForSubTask(const misc_htm_manager_t manager, misc_size_t taskId, misc_size_t percentComplete);
int misc_htm_runSubTask(misc_htm_manager_t restrict manager, misc_size_t taskId, misc_htm_subTaskFunction_t subTask,
                        void** restrict data, misc_size_t numPieces);

misc_size_t misc_htm_getNumThreadsForTopLevelTask(const misc_htm_manager_t threadManager, misc_size_t taskId);
void misc_htm_getNumPiecesForSubTask(const misc_htm_manager_t restrict threadManager, misc_size_t taskId,
                                     misc_size_t numElements, misc_size_t minNumElementsPerPiece,
                                     misc_size_t* restrict numPiecesPtr, misc_size_t* restrict numElementsPerPiecePtr,
                                     misc_size_t* restrict offByOneIndexPtr);

void misc_htm_printf(misc_htm_manager_t manager, const char* format, ...);

// blocking thread manager
#define EXT_BTM_INVALID_THREAD_ID ((misc_size_t) -1)
struct _misc_btm_manager_t;
typedef struct _misc_btm_manager_t* misc_btm_manager_t;
typedef void (*misc_btm_taskFunction_t)(void*);

int misc_btm_create(misc_btm_manager_t* manager, misc_size_t numThreads);
int misc_btm_destroy(misc_btm_manager_t manager);

int misc_btm_runTasks(misc_btm_manager_t restrict manager, misc_btm_taskFunction_t task,
                      void** restrict data, misc_size_t numTasks);
misc_size_t misc_btm_getNumThreads(const misc_btm_manager_t manager);
void misc_btm_getNumThreadsForJob(const misc_btm_manager_t restrict threadManager, misc_size_t numElements, misc_size_t minNumElementsPerThread,
                                  misc_size_t* restrict numThreadsPtr, misc_size_t* restrict numElementsPerThreadPtr,
                                  misc_size_t* restrict offByOneIndexPtr);

bool misc_btm_isNull(misc_btm_manager_t manager);

int misc_btm_runTaskInParentThread(misc_btm_manager_t restrict manager, misc_size_t threadId, misc_btm_taskFunction_t task,
                                   void* restrict data);


#ifdef __cplusplus
}
#endif

#endif // MISC_THREAD_H

