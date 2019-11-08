#ifndef STIMER
#define STIMER

#include <cassert>
#include <time.h>

// We sometimes make arrays of STimers.  Let's avoid false sharing!
#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 128
#endif

class alignas(CACHE_LINE_SIZE) STimer {
public:
    STimer();
    ~STimer();
    void Start(void);
    void Stop(void);
    double Elapsed(void);
    void Clear(void);
    void increment(struct timespec dt);
    struct timespec get_timer(void);
    int timeron;

    struct timespec timer;
    
private:
    struct timespec tstart;
};

struct timespec scale_timer(double s, struct timespec t);
void timespecclear(struct timespec *t);
void timespecadd(struct timespec *a, struct timespec *b, struct timespec *res);
void timespecsub(struct timespec *a, struct timespec *b, struct timespec *res);
#endif
