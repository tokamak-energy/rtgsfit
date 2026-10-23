#define _POSIX_C_SOURCE 199309L

#include "timing.h"

#include <time.h>

static uint64_t timer_total_ns[RTGSFIT_TIMER_COUNT];

#ifdef RTGSFIT_TIMING

uint64_t timing_now_ns(void)
{
    struct timespec now;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &now) != 0) {
        return 0;
    }
    return (uint64_t)now.tv_sec * 1000000000u + (uint64_t)now.tv_nsec;
}

uint64_t timing_accumulate(enum rtgsfit_timer timer, uint64_t start_ns)
{
    const uint64_t now_ns = timing_now_ns();
    timer_total_ns[timer] += now_ns - start_ns;
    return now_ns;
}

#endif /* RTGSFIT_TIMING */

void rtgsfit_timing_reset(void)
{
    for (int32_t i = 0; i < RTGSFIT_TIMER_COUNT; i++) {
        timer_total_ns[i] = 0;
    }
}

int32_t rtgsfit_timing_get(double time_us[], int32_t n_time_us)
{
    if (n_time_us != RTGSFIT_TIMER_COUNT) {
        return RTGSFIT_ERR_ARGUMENT_SIZE;
    }
    for (int32_t i = 0; i < RTGSFIT_TIMER_COUNT; i++) {
        time_us[i] = (double)timer_total_ns[i] * 1e-3;
    }
    return RTGSFIT_OK;
}
