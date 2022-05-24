#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "decs.h"

int get_wctime(double * wclk) {
    struct timeval mytime;
    if(gettimeofday(&mytime, NULL)) {
        return ERROR_INTERNAL;
    }
    *wclk = (double)mytime.tv_sec + (double)mytime.tv_usec*1e-6;
    return ERROR_NONE;
}

int timer_init(char const * name, Timer ** timer) {
    Timer * out = malloc(sizeof(*out));

    out->name = malloc(strlen(name) + 1);
    strcpy((char *)out->name, name);

    out->running = false;
    out->tick = 0.0;
    out->tock = 0.0;
    out->ncycles = 0;
    out->total = 0.0;

    *timer = out;

    return ERROR_NONE;
}

int timer_start(Timer * timer) {
    if(false == timer->running) {
        int ierr = get_wctime(&timer->tick); CHKERRQ(ierr);
        timer->running = true;
    }
    return ERROR_NONE;
}

int timer_stop(Timer * timer) {
    int ierr = get_wctime(&timer->tock); CHKERRQ(ierr);
    timer->running = false;
    timer->total += timer->tock - timer->tick;
    ++timer->ncycles;
    return ERROR_NONE;
}

int timer_reset(Timer * timer) {
    if(timer->running) {
        return ERROR_INTERNAL;
    }
    timer->total = 0;
    timer->tick = 0;
    timer->tock = 0;
    timer->ncycles = 0;
    return ERROR_NONE;
}

int timer_get_ncycles(Timer const * timer, long unsigned int * ncycles) {
    *ncycles = timer->ncycles;
    return ERROR_NONE;
}

int timer_get_total(Timer const * timer, double * total) {
    *total = timer->total;
    return ERROR_NONE;
}

int timer_get_name(Timer const * timer, char const ** name) {
    *name = timer->name;
    return ERROR_NONE;
}

int timer_free(Timer * timer) {
    free((char *)timer->name);
    free(timer);
    return ERROR_NONE;
}
