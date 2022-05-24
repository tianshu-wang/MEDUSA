//! \file timer.h Simple timer class

#ifndef TIMER_H
#define TIMER_H

#include <stdbool.h>

int get_wctime(double * wckl);

//! Simple timer class
typedef struct {
    //! True if the clock is running
    bool running;
    //! Start timestamp
    double tick;
    //! Stop timestamp
    double tock;
    //! Number of cycles
    long unsigned int ncycles;
    //! Total elapsed time
    double total;
    //! Timer name
    char const * name;
} Timer;

//! Initialize a timer
int timer_init(
        //! [in] Name of the timer
        char const * name,
        //! [out] Newly allocated timer
        Timer ** timer);

//! Start a timer
int timer_start(
        //! [in,out] Timer to start
        Timer * timer);

//! Stop the timer
int timer_stop(
        //! [in,out] Timer to start
        Timer * timer);

//! Reset the timer
int timer_reset(
        //! [in,out] Timer to reset
        Timer * timer);

//! Total number of cycles
int timer_get_ncycles(
        //! [in] Timer object
        Timer const * timer,
        //! [out] Total elapsed time
        long unsigned int * ncycles);

//! Get the total elapsed time
int timer_get_total(
        //! [in] Timer object
        Timer const * timer,
        //! [out] Total elapsed time
        double * total);

//! Get the timer name
int timer_get_name(
        //! [in] Timer object
        Timer const * timer,
        //! [out] Timer name
        char const ** name);

//! De-allocate the timer
int timer_free(
        //! [in,out] Timer to de-allocate
        Timer * timer);

#endif

/* vim: set ft=c.doxygen : */
