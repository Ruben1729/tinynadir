/* SPDX-License-Identifier: MIT */
#ifndef TINYPKEPLER_H
#define TINYPKEPLER_H
#include "tinyelements/tinyelements.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYPKEPLER_OK,TINYPKEPLER_INVALID,TINYPKEPLER_NO_CONVERGENCE } tinypkepler_status;
typedef struct {tinyelements_classical mean_elements;tinyelements_state state;double raan_rate_rad_s,arg_periapsis_rate_rad_s,mean_motion_rad_s;} tinypkepler_result;
/* Vallado-style pkepler J2 secular-precession profile, ndot=nddot=0.
   Elliptic model-specific mean elements (p>0, 0<=e<1), not SGP4/TLE elements.
   p,a,e,i remain constant; node/periapsis use first-order J2 rates and
   mean anomaly advances at Keplerian sqrt(mu/a^3), without a J2 mean-rate correction.
   SI units and radians. No J4, drag-rate, short-period reconstruction, or tides.
   |n*elapsed_s|<=1e8. Failures preserve output. */
tinypkepler_status tinypkepler_propagate(double mu,double equatorial_radius_m,double j2,
 const tinyelements_classical *initial_mean_elements,double elapsed_s,tinypkepler_result *output);
#ifdef __cplusplus
}
#endif
#endif
