#ifndef EXTERNAL_THREAD_H
#define EXTERNAL_THREAD_H

#include "stddef.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _ext_mt_manager_t;
typedef struct _ext_mt_manager_t* ext_mt_manager_t;


int ext_mt_create(ext_mt_manager_t* manager, ext_size_t numThreads);
int ext_mt_destroy(ext_mt_manager_t manager);

ext_size_t ext_mt_getNumThreads(const ext_mt_manager_t manager);
void ext_mt_getNumThreadsForJob(const ext_mt_manager_t restrict threadManager, ext_size_t numElements,
                                ext_size_t minNumElementsPerThread, ext_size_t* restrict numThreadsPtr, ext_size_t* restrict numElementsPerThreadPtr);


typedef void (*ext_mt_taskFunction_t)(void*);
  
int ext_mt_runTasks(ext_mt_manager_t restrict manager, ext_mt_taskFunction_t task,
                    void** restrict data, ext_size_t numTasks);


#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_THREAD_H
