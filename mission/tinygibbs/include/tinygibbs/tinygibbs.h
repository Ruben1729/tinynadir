/* SPDX-License-Identifier: MIT */
#ifndef TINYGIBBS_H
#define TINYGIBBS_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYGIBBS_OK, TINYGIBBS_INVALID, TINYGIBBS_GEOMETRY } tinygibbs_status;
typedef struct {double velocity_m_s[3],angle12_rad,angle23_rad,coplanarity_rad;} tinygibbs_result;
/* Three consecutive position observations in one inertial frame, metres.
   mu m^3/s^2; caller chooses accepted coplanarity error in [0,pi/2].
   Velocity belongs to the middle position. Nearly collinear/degenerate geometry
   is rejected; close arcs are better served by Herrick-Gibbs. Failure preserves output. */
tinygibbs_status tinygibbs_estimate(double mu,const double r1_m[3],const double r2_m[3],
                                  const double r3_m[3],double max_coplanarity_rad,tinygibbs_result *output);
#ifdef __cplusplus
}
#endif
#endif
