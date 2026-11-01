/* SPDX-License-Identifier: MIT */
#ifndef TINYSGP4FIT_H
#define TINYSGP4FIT_H
#include <stddef.h>
#include "tinytle/tinytle.h"
#include "tinybatchls/tinybatchls.h"
#ifdef __cplusplus
extern "C" {
#endif
#define TINYSGP4FIT_MAX_EPOCHS 10
typedef enum {TINYSGP4FIT_OK,TINYSGP4FIT_INVALID,TINYSGP4FIT_MODEL_ERROR,TINYSGP4FIT_UNOBSERVABLE,TINYSGP4FIT_NO_CONVERGENCE} tinysgp4fit_status;
typedef struct {
    double minutes_since_epoch;
    double position_m[3],velocity_m_s[3];
    double position_sigma_m[3],velocity_sigma_m_s[3];
} tinysgp4fit_observation;
typedef struct {size_t max_iterations;int fit_bstar,use_velocity;double initial_damping;} tinysgp4fit_options;
typedef struct {
    tinytle_elements elements;
    /* Packed parameter_count^2 row-major covariance. Parameter order/units:
       inclination[deg], RAAN[deg], eccentricity, argument of perigee[deg],
       mean anomaly[deg], mean motion[rev/day], and optional B*[1/Earth radius]. */
    double covariance[49],weighted_sum_squares,normalized_rms;
    size_t parameter_count,scalar_observation_count,iterations,model_evaluations;
} tinysgp4fit_result;
typedef struct {
    tinybatchls_workspace batch;
    double scalar_observations[64],sigma[64];
    tgp4_satellite satellite;
    tgp4_state state;
} tinysgp4fit_workspace;
tinysgp4fit_status tinysgp4fit_default_options(tinysgp4fit_options *output);
/* Local weighted differential correction of conventional SGP4 GP elements.
   Epoch, catalog and ephemeris metadata remain fixed. Observations must be
   same-frame TEME positions/coordinate velocities, relative to the seed epoch.
   This profile fits 6 parameters, or 7 including B*. It does not fit epoch,
   ndot/nddot, light-time, station biases, maneuvers or XP elements.
   The classical parameterization excludes e<=1e-12 or e>=1-1e-12 and inclinations
   within1e-6 deg of0/180. Exactly circular/equatorial element fits are unobservable.
   Supply2..10 epochs (>=3 for position-only), times within +/-10080 minutes,
   positive independent sigmas for included components, and a nearby seed.
   Scratch must not overlap inputs/output and is modified on failure.
   Result is preserved on failure. No heap, file I/O, or mutable globals. */
tinysgp4fit_status tinysgp4fit_fit(const tinytle_elements *seed,
 const tinysgp4fit_observation *observations,size_t epoch_count,
 const tinysgp4fit_options *options,tinysgp4fit_workspace *workspace,tinysgp4fit_result *output);
#ifdef __cplusplus
}
#endif
#endif
