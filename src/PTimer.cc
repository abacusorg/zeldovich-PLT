/* PTimer.cc - Template implementation (included by PTimer.h, do not include directly) */
#ifndef PTIMER_CC_IMPL
#define PTIMER_CC_IMPL

/* Parallel Timer - Template implementation
 * 
 * PTimer is a high-resolution, low-overhead timer for performance-intensive
 * parallel regions of the code.  We try to detect the host architecture and
 * use the fastest available timer for that platform.  Usually this means
 * that the resulting times are not wall-clock times, but instead something
 * like CPU cycles.
 *
 * If you need wall-clock time, the PTimerWall class is provided.  One could
 * use this for timing multi-threaded IO, for example. Of course, one can always
 * use an array of STimer for that kind of task, but then you have to manage
 * the array memory allocation manually.
 */

#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include "PTimer.h"

template <class HardwareTimer>
PTimerBase<HardwareTimer>::PTimerBase(void) : PTimerBase(omp_get_max_threads()) { }

template <class HardwareTimer>
PTimerBase<HardwareTimer>::PTimerBase(int nthreads) { 
    nprocs = nthreads;

    assert(posix_memalign((void **) &timer, CACHE_LINE_SIZE, sizeof(*timer)*nprocs) == 0);

    for(int g=0; g<nprocs; g++)
        timer[g].on = 0;
    Clear();
}

template <class HardwareTimer>
PTimerBase<HardwareTimer>::~PTimerBase() {
    free(timer);
}

template <class HardwareTimer>
void PTimerBase<HardwareTimer>::Clear(void) {
    for(int g=0; g<nprocs; g++)
        assert(!timer[g].on);

    for(int g=0; g<nprocs; g++) {
        timer[g].Clear();
    }
}

template <class HardwareTimer>
inline void PTimerBase<HardwareTimer>::Start(void){
    Start(omp_get_thread_num());
}

template <class HardwareTimer>
inline void PTimerBase<HardwareTimer>::Start(int thread_num) {
    int g = thread_num;
    
    assert(g < nprocs);  // If this fails omp_get_max_threads() may not be returning the global max # of threads
    assert(!timer[g].on);

    timer[g].Start();
    timer[g].on = 1;

#ifdef PTIMER_PPC
    tbr_count = __ppc_get_timebase();
#endif
}

template <class HardwareTimer>
inline void PTimerBase<HardwareTimer>::Stop(void){
    Stop(omp_get_thread_num());
}

template <class HardwareTimer>
inline void PTimerBase<HardwareTimer>::Stop(int thread_num) {
    int g = thread_num;
    assert(timer[g].on);

    timer[g].Stop();
    timer[g].on = 0;
}

template <class HardwareTimer>
double PTimerBase<HardwareTimer>::Elapsed(void) {
    for(int g=0; g<nprocs; g++) 
        assert(!timer[g].on);  
    
    double sum = 0;

    for(int g=0; g<nprocs; g++)
        sum += timer[g].Elapsed();

    return sum;
}


// Now provide template specializations for a dummy timer
template <>
inline PTimerBase<DummyHardwareTimer>::PTimerBase(int nthreads [[maybe_unused]]) { }
template <>
inline PTimerBase<DummyHardwareTimer>::~PTimerBase() { }
template <>
inline void PTimerBase<DummyHardwareTimer>::Start(void){ }
template <>
inline void PTimerBase<DummyHardwareTimer>::Start(int thread_num [[maybe_unused]]) { }
template <>
inline void PTimerBase<DummyHardwareTimer>::Stop(void){ }
template <>
inline void PTimerBase<DummyHardwareTimer>::Stop(int thread_num [[maybe_unused]]) { }
template <>
inline double PTimerBase<DummyHardwareTimer>::Elapsed(void){ return 1e-12; }
    // We return a small number so that we don't divide by zero

#endif /* PTIMER_CC_IMPL */
