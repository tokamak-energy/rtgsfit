#include "measurement_processing.h"
#include "test_utilities.h"

static void test_prepare(void)
{
    const double meas[] = {2.0, 4.0, 8.0};
    const double coil_current[] = {3.0};
    const double replacement[] = {1.0, 0.0, 0.0, 0.0, 0.5, 0.5};
    const double g_coil[] = {2.0, -1.0};
    const double weight[] = {0.5, 2.0, 7.0};
    double measured[2], coil[2], fit_target[3];

    CHECK(measurement_processing_prepare(meas, 3, coil_current, 1, replacement, 6, g_coil, 2,
                                         weight, 3, measured, 2, coil, 2, fit_target,
                                         3) == RTGSFIT_OK);
    CHECK_CLOSE(measured[0], 2.0, 0.0);
    CHECK_CLOSE(measured[1], 6.0, 0.0);
    CHECK_CLOSE(coil[0], 6.0, 0.0);
    CHECK_CLOSE(coil[1], -3.0, 0.0);
    CHECK_CLOSE(fit_target[0], -2.0, 0.0);
    CHECK_CLOSE(fit_target[1], 18.0, 0.0);
    CHECK_CLOSE(fit_target[2], 0.0, 0.0);
}

static void test_prepare_rejects_invalid_buffers(void)
{
    const double meas[3] = {0.0};
    const double coil_current[1] = {0.0};
    const double replacement[6] = {0.0};
    const double g_coil[2] = {0.0};
    const double weight[3] = {0.0};
    double measured[2], coil[2], fit_target[3];

    for (int32_t i = 0; i < 8; i++)
    {
        for (int32_t delta = -1; delta <= 1; delta += 2)
        {
            int32_t sizes[] = {3, 1, 6, 2, 3, 2, 2, 3};
            sizes[i] += delta;
            measured[0] = coil[0] = fit_target[0] = 123.0;
            CHECK(measurement_processing_prepare(
                      meas, sizes[0], coil_current, sizes[1], replacement, sizes[2], g_coil,
                      sizes[3], weight, sizes[4], measured, sizes[5], coil, sizes[6], fit_target,
                      sizes[7]) == RTGSFIT_ERR_ARGUMENT_SIZE);
            CHECK_CLOSE(measured[0], 123.0, 0.0);
            CHECK_CLOSE(coil[0], 123.0, 0.0);
            CHECK_CLOSE(fit_target[0], 123.0, 0.0);
        }
    }
    CHECK(measurement_processing_prepare(NULL, 3, coil_current, 1, replacement, 6, g_coil, 2,
                                         weight, 3, measured, 2, coil, 2, fit_target,
                                         3) == RTGSFIT_ERR_ARGUMENT_SIZE);
}

int main(void)
{
    RUN_TEST(test_prepare);
    RUN_TEST(test_prepare_rejects_invalid_buffers);
    return test_failures == 0 ? 0 : 1;
}
