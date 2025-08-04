#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

namespace HexHex::Threading
{

static inline void setNumThreads(const size_t num_threads)
{
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
}

static inline size_t getNumThreads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

static inline size_t getThreadID()
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

}

