#ifndef MISC_PTHREAD_H
#define MISC_PTHREAD_H

#include <pthread.h>

typedef pthread_t Thread;
typedef pthread_mutex_t Mutex;
typedef pthread_cond_t Condition;

#define createThread(_THREAD_, _FUNC_, _DATA_) pthread_create(_THREAD_, NULL, _FUNC_, _DATA_)
#define joinThread(_X_) pthread_join(_X_, NULL)

#define initializeMutex(_X_) pthread_mutex_init(&(_X_), NULL)
#define destroyMutex(_X_) pthread_mutex_destroy(&(_X_))
#define lockMutex(_X_) pthread_mutex_lock(&(_X_))
#define unlockMutex(_X_) pthread_mutex_unlock(&(_X_))

#define initializeCondition(_X_) pthread_cond_init(&(_X_), NULL)
#define destroyCondition(_X_) pthread_cond_destroy(&(_X_))
#define waitOnCondition(_COND_, _MUTEX_) pthread_cond_wait(&(_COND_), &(_MUTEX_))
#define signalCondition(_X_) pthread_cond_signal(&(_X_))
#define waitOnConditionForTime(_COND_, _MUTEX_, _TIME_) pthread_cond_timedwait(&(_COND_), &(_MUTEX_), &(_TIME_))

#endif

