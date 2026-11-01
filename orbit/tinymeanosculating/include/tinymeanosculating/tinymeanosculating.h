/* SPDX-License-Identifier: MIT */
#ifndef TINYMEANOSCULATING_H
#define TINYMEANOSCULATING_H
#include <stddef.h>
#include "tinycowell/tinycowell.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYMEANOSCULATING_OK=0,TINYMEANOSCULATING_INVALID,TINYMEANOSCULATING_RANGE,TINYMEANOSCULATING_BUDGET,TINYMEANOSCULATING_INTEGRATION,TINYMEANOSCULATING_NO_CONVERGENCE} tinymeanosculating_status;
typedef struct {double mu_m3_s2,equatorial_radius_m,j2;} tinymeanosculating_model;
/* Numerically averaged prograde equinoctial mean coordinates:
   a, f=e*cos(Omega+omega), g=e*sin(Omega+omega), h=tan(i/2)*cos(Omega),
   k=tan(i/2)*sin(Omega), lambda=Omega+omega+M at the center epoch.
   This module's centered-one-orbit time-average definition is NOT Brouwer,
   DSST, SGP4 or TLE mean elements. See README for the precise filter. */
typedef struct {double semimajor_axis_m,f,g,h,k,mean_longitude_rad;} tinymeanosculating_mean;
typedef struct {
 unsigned samples; /* even,16..512 trapezoid intervals over one centered orbit */
 double max_step_s; /* positive RK4 step bound */
 size_t max_steps; /* total integration-step budget,1..1000000, entire call */
 unsigned max_iterations; /* inverse averaging passes,1..16 */
 double inverse_tolerance; /* scaled a residual and absolute other fields,1e-12..1e-4 */
} tinymeanosculating_options;
typedef struct {double integration[48];} tinymeanosculating_workspace;
typedef struct {size_t integration_steps;unsigned averaging_passes;double residual;} tinymeanosculating_stats;
/* Model profile: central+J2 gravity, inertial fixed +z pole, |J2|<=.005;
   elliptic e<=.3, inclination<=170deg and osculating pericenter>body radius.
   SI/radians. Caller owns workspace, disjoint from inputs/outputs. No heap.
   Output/stats preserved on failure; workspace unspecified. */
tinymeanosculating_status tinymeanosculating_from_osculating(const tinymeanosculating_model *model,const tinycowell_state *state,const tinymeanosculating_options *options,tinymeanosculating_workspace *workspace,tinymeanosculating_mean *out,tinymeanosculating_stats *stats);
/* Bounded fixed-point inverse of the SAME numerical filter. Initial candidate
   treats mean coordinates as osculating, then corrects averaged residuals.
   Can fail to converge; never silently returns an uncorrected mean state. */
tinymeanosculating_status tinymeanosculating_to_osculating(const tinymeanosculating_model *model,const tinymeanosculating_mean *mean,const tinymeanosculating_options *options,tinymeanosculating_workspace *workspace,tinycowell_state *out,tinymeanosculating_stats *stats);
#ifdef __cplusplus
}
#endif
#endif
