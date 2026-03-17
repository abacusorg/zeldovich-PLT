/* PTimer.cc
 * 
 * Parallel Timer
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
 *
 * Define the PTIMER_DUMMY macro before #including this file to turn off
 * the timers throughout the code (PTimerWall will continue to function).
 */

#ifndef PTIMER
#define PTIMER

#include <cassert>
#include <cstdint>
#include <time.h>

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 128
#endif

typedef uint64_t uint64;

template <class HardwareTimer>
class PTimerBase {
public:
    PTimerBase(void);
    PTimerBase(int nthreads);
    ~PTimerBase();

    void Start(void);
    void Stop(void);
    void Start(int thread_num);
    void Stop(int thread_num);
    double Elapsed(void);
    void Clear(void);

    int nprocs;
    HardwareTimer *timer;
};


// A portable implementation of HardwareTimer
// This also uses seconds/nanoseconds, so we can use it for wall-clock time
struct alignas(CACHE_LINE_SIZE) timespec_timer {
    struct timespec tstart;
    struct timespec tot;
    int on;

    void Start(void){
        assert( clock_gettime( CLOCK_MONOTONIC, &tstart) == 0 );
    }

    void Stop(void){
        struct timespec dt;
        struct timespec *t = &tot;
        assert( clock_gettime( CLOCK_MONOTONIC, &dt) == 0 );
        t->tv_nsec += dt.tv_nsec - tstart.tv_nsec;
        t->tv_sec += dt.tv_sec - tstart.tv_sec;
        if (t->tv_nsec<0) { t->tv_nsec += 1000000000; t->tv_sec--; }
        else if (t->tv_nsec>1000000000) { t->tv_nsec -= 1000000000; t->tv_sec++; }
    }

    double Elapsed(void){
        return tot.tv_sec + 1e-9*tot.tv_nsec;
    }

    void Clear(void){
        tot.tv_sec = 0;
        tot.tv_nsec = 0;
    }
};

struct DummyHardwareTimer { };


#ifndef PTIMER_DUMMY
// PTimer is expensive enough (estimate is about 35 ns per call on 'ted')
// to notice, because it's used only in parallel loops where we're often
// working on relatively small chunks.  So we provide a way to no-op these.

// Use the fastest timing function available on the platform
#ifdef _ARCH_PPC

#include <sys/platform/ppc.h>
struct alignas(CACHE_LINE_SIZE) tbr_timer {
    uint64 tbr_start;
    uint64 tbr_tot;
    int on;

    void Start(void){
        tbr_start = __ppc_get_timebase();
    }

    void Stop(void){
        tbr_tot += __ppc_get_timebase() - tbr_start;
    }

    uint64 Elapsed(void){
        return tbr_tot;
    }

    void Clear(void){
        tbr_start = 0;
        tbr_tot = 0;
    }
};

typedef PTimerBase<struct tbr_timer> PTimer;

#elif defined(__x86_64)

#include <x86intrin.h>
struct alignas(CACHE_LINE_SIZE) tsc_timer {
    uint64 tsc_start;
    uint64 tsc_tot;
    int on;

    void Start(void){
        tsc_start = __rdtsc();
    }

    void Stop(void){
        tsc_tot += __rdtsc() - tsc_start;
    }

    uint64 Elapsed(void){
        return tsc_tot;
    }

    void Clear(void){
        tsc_start = 0;
        tsc_tot = 0;
    }
};

typedef PTimerBase<struct tsc_timer> PTimer;

#else

typedef PTimerBase<struct timespec_timer> PTimer;

#endif

#else // PTIMER_DUMMY

typedef PTimerBase<DummyHardwareTimer> PTimer;

#endif // PTIMER_DUMMY

// We probably want the PTimerWall timers to keep working even with PTIMER_DUMMY
typedef PTimerBase<struct timespec_timer> PTimerWall;

// Include template implementation (PTimer is header-only style for templates)
#include "PTimer.cc"

#endif
