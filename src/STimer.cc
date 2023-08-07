#include <cstdio>
#include <cstdlib>
#include "STimer.h"

STimer::STimer(void) { 
    timespecclear(&timer);
    timeron = 0;
}

STimer::~STimer() {
}

struct timespec STimer::get_timer(void) {
    return timer;
}

void STimer::increment(struct timespec dt) {
    timespecadd(&dt, &timer, &timer);
}

void STimer::Start() {
    assert(!timeron); 
    assert( clock_gettime( CLOCK_MONOTONIC, &tstart) == 0 );
    timeron = 1;
}

void STimer::Stop(void) {
    assert( timeron );
    struct timespec dt, tend;
    assert( clock_gettime( CLOCK_MONOTONIC, &tend) == 0 );
    timespecsub(&tend, &tstart, &dt);
    timespecadd(&dt, &timer, &timer);
    timeron = 0;
}

void STimer::Clear(void) {
    assert(!timeron); 
    timespecclear(&timer);
}

double STimer::Elapsed(void) {
    return  timer.tv_sec + 1e-9*timer.tv_nsec;
}

/* Define timespec addition, subtraction, scaling, resetting.
 * We assume that timespec are well behaved, with tv_nsec values
 * between 0 and 1 billion.
 */

#define NSEC_PER_SEC 1000000000

struct timespec scale_timer(double s, struct timespec t) {
    int64_t ns = t.tv_sec*NSEC_PER_SEC + t.tv_nsec;
    ns *= s;

    struct timespec tp;
    tp.tv_sec = ns / NSEC_PER_SEC;
    tp.tv_nsec = ns - tp.tv_sec*NSEC_PER_SEC;
    return tp;
}

void timespecclear(struct timespec *t){
    t->tv_sec = 0;
    t->tv_nsec = 0;
}

void timespecadd(struct timespec *a, struct timespec *b, struct timespec *res){
    res->tv_sec = a->tv_sec + b->tv_sec;
    res->tv_nsec = a->tv_nsec + b->tv_nsec;
    if(res->tv_nsec >= NSEC_PER_SEC){
        res->tv_sec += 1;
        res->tv_nsec -= NSEC_PER_SEC;
    }
}

void timespecsub(struct timespec *a, struct timespec *b, struct timespec *res){
    res->tv_sec = a->tv_sec - b->tv_sec;
    res->tv_nsec = a->tv_nsec - b->tv_nsec;

    // Could use the following if one is paranoid about unsigned time_t
    //if(a->tv_sec < b->tv_sec && res->tv_sec > 0){
    //    assert(0 && "struct timespec.tv_sec must be signed type!");
    //}

    // Allow negative delta_t
    // This could result in invalid timespecs, though!
    if(res->tv_nsec < 0 && res->tv_sec > 0){
        res->tv_sec -= 1;
        res->tv_nsec += NSEC_PER_SEC;
    }
    else if(res->tv_nsec > 0 && res->tv_sec < 0){
        res->tv_sec += 1;
        res->tv_nsec -= NSEC_PER_SEC;
    }
}
