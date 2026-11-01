/* SPDX-License-Identifier: MIT */
#ifndef TINYHERRICKGIBBS_H
#define TINYHERRICKGIBBS_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYHERRICKGIBBS_OK, TINYHERRICKGIBBS_INVALID, TINYHERRICKGIBBS_GEOMETRY } tinyherrickgibbs_status;
typedef struct {double velocity_m_s[3],angle12_rad,angle23_rad,coplanarity_rad;} tinyherrickgibbs_result;
/* Short-arc velocity at middle observation. Positions metres in one inertial
   frame; mu m^3/s^2. dt_before=t2-t1>0, dt_after=t3-t2>0, seconds.
   Supply differences rather than subtracting large Julian dates internally.
   Caller selects max adjacent angle and max coplanarity error, radians.
   Each tolerance in (0,pi/2]; failures preserve output. */
tinyherrickgibbs_status tinyherrickgibbs_estimate(double mu,const double r1_m[3],
 const double r2_m[3],const double r3_m[3],double dt_before_s,double dt_after_s,
 double max_angle_rad,double max_coplanarity_rad,tinyherrickgibbs_result *output);
#ifdef __cplusplus
}
#endif
#endif
