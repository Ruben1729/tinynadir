#ifndef TINYSLEW_H
#define TINYSLEW_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYSLEW_OK=0,TINYSLEW_INVALID_ARGUMENT,TINYSLEW_DEGENERATE,TINYSLEW_ARITHMETIC_ERROR} tinyslew_status;
typedef struct {double initial[4],target[4],axis[3],angle,acceleration,peak_rate,accel_time,cruise_time,duration;} tinyslew_plan;
/* Shortest geodesic, rest-to-rest, time-minimal triangular/trapezoidal angular
 * speed profile. q is scalar-first Hamilton reference->body; maximum_rate in
 * rad/s and maximum_acceleration in rad/s² are positive. Not jerk-limited.
 * A plan is immutable after successful creation. Failed calls preserve outputs. */
tinyslew_status tinyslew_initialize(tinyslew_plan *plan,const double initial[4],const double target[4],double maximum_rate,double maximum_acceleration);
/* elapsed>=0 seconds; after completion returns target/zero rate/acceleration.
 * Outputs are q and physical desired angular rate/acceleration in current
 * desired body axes. Output arrays must be distinct. */
tinyslew_status tinyslew_sample(const tinyslew_plan *plan,double elapsed,double q[4],double rate[3],double acceleration[3]);
#ifdef __cplusplus
}
#endif
#endif
