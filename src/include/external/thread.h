#ifndef EXTERNAL_THREAD_H
#define EXTERNAL_THREAD_H

#include "stddef.h"
#include <stdbool.h>
#define __USE_POSIX199309 1
#include <time.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _ext_mt_manager_t;
typedef struct _ext_mt_manager_t* ext_mt_manager_t;


int ext_mt_create(ext_mt_manager_t* manager, ext_size_t numThreads);
int ext_mt_destroy(ext_mt_manager_t manager);

ext_size_t ext_mt_getNumThreads(const ext_mt_manager_t manager);
// assign numElementsPerThread for the first offByOneIndex threads, and numElementsPerThread - 1 for those after
void ext_mt_getNumThreadsForJob(const ext_mt_manager_t restrict threadManager, ext_size_t numElements, ext_size_t minNumElementsPerThread,
                                ext_size_t* restrict numThreadsPtr, ext_size_t* restrict numElementsPerThreadPtr, ext_size_t* restrict offByOneIndexPtr);


typedef void (*ext_mt_taskFunction_t)(void*);

typedef void (*ext_mt_infoFunction_t)(void** threadData, ext_size_t numThreads);

int ext_mt_runTasks(ext_mt_manager_t restrict manager, ext_mt_taskFunction_t task,
                    void** restrict data, ext_size_t numTasks);

int ext_mt_runTasksWithInfo(ext_mt_manager_t restrict manager, ext_mt_taskFunction_t task,
                    void** restrict data, ext_size_t numTasks, time_t sleepSeconds, ext_mt_infoFunction_t info);


// hierarchical thread manager
#define EXT_HTM_INVALID_TASK_ID ((ext_size_t) -1)
struct _ext_htm_manager_t;
typedef struct _ext_htm_manager_t* ext_htm_manager_t;


typedef void (*ext_htm_topLevelTaskFunction_t)(ext_size_t taskId, void*);
typedef void (*ext_htm_subTaskFunction_t)(void*);

int ext_htm_create(ext_htm_manager_t* manager, ext_size_t numThreads);
int ext_htm_destroy(ext_htm_manager_t manager);

int ext_htm_runTopLevelTasks(ext_htm_manager_t restrict manager, ext_htm_topLevelTaskFunction_t task,
                             void** restrict data, ext_size_t numTasks);
int ext_htm_runTopLevelTasksWithOutput(ext_htm_manager_t restrict manager, ext_htm_topLevelTaskFunction_t task,
                                       void** restrict data, ext_size_t numTasks, const struct timespec* restrict outputDelay);

ext_size_t ext_htm_reserveThreadsForSubTask(const ext_htm_manager_t manager, ext_size_t taskId, ext_size_t percentComplete);
int ext_htm_runSubTask(ext_htm_manager_t restrict manager, ext_size_t taskId, ext_htm_subTaskFunction_t subTask,
                       void** restrict data, ext_size_t numPieces);

ext_size_t ext_htm_getNumThreadsForTopLevelTask(const ext_htm_manager_t threadManager, ext_size_t taskId);
void ext_htm_getNumPiecesForSubTask(const ext_htm_manager_t restrict threadManager, ext_size_t taskId,
                                    ext_size_t numElements, ext_size_t minNumElementsPerPiece,
                                    ext_size_t* restrict numPiecesPtr, ext_size_t* restrict numElementsPerPiecePtr, ext_size_t* restrict offByOneIndexPtr);

void ext_htm_printf(ext_htm_manager_t manager, const char* format, ...);

// blocking thread manager
#define EXT_BTM_INVALID_THREAD_ID ((ext_size_t) -1)
struct _ext_btm_manager_t;
typedef struct _ext_btm_manager_t* ext_btm_manager_t;
typedef void (*ext_btm_taskFunction_t)(void*);

int ext_btm_create(ext_btm_manager_t* manager, ext_size_t numThreads);
int ext_btm_destroy(ext_btm_manager_t manager);

int ext_btm_runTasks(ext_btm_manager_t restrict manager, ext_btm_taskFunction_t task,
                     void** restrict data, ext_size_t numTasks);
ext_size_t ext_btm_getNumThreads(const ext_btm_manager_t manager);
void ext_btm_getNumThreadsForJob(const ext_btm_manager_t restrict threadManager, ext_size_t numElements, ext_size_t minNumElementsPerThread,
                                 ext_size_t* restrict numThreadsPtr, ext_size_t* restrict numElementsPerThreadPtr, ext_size_t* restrict offByOneIndexPtr);

bool ext_btm_isNull(ext_btm_manager_t manager);

int ext_btm_runTaskInParentThread(ext_btm_manager_t restrict manager, ext_size_t threadId, ext_btm_taskFunction_t task, void* restrict data);


#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_THREAD_H
