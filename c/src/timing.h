/*
 * Per-section CPU timing of `rtgsfit()`, enabled by building with `TIMING=1`
 * (which defines RTGSFIT_TIMING). When it is disabled, the calls compile to nothing.
 */
#ifndef RTGSFIT_TIMING_H
#define RTGSFIT_TIMING_H

#include <stdint.h>

#include "rtgsfit.h"

#ifdef RTGSFIT_TIMING

/* Thread CPU time in nanoseconds. */
uint64_t timing_now_ns(void);

/* Add the time since `start_ns` to `timer`, and return the current time. */
uint64_t timing_accumulate(enum rtgsfit_timer timer, uint64_t start_ns);

#else

static inline uint64_t timing_now_ns(void)
{
    return 0;
}

static inline uint64_t timing_accumulate(enum rtgsfit_timer timer, uint64_t start_ns)
{
    (void)timer;
    (void)start_ns;
    return 0;
}

#endif /* RTGSFIT_TIMING */

#endif /* RTGSFIT_TIMING_H */
