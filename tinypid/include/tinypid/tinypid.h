#ifndef TINYPID_H
#define TINYPID_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYPID_OK=0,TINYPID_INVALID_ARGUMENT,TINYPID_DEGENERATE,TINYPID_ARITHMETIC_ERROR} tinypid_status;
typedef struct {double integral[3];} tinypid_state;
/* Zero-initialize integral. Uses tinypd's reference->body convention and rate
 * error (measured-desired in current body axes). Integrates 2*sign(qerr.w)*qerr.v
 * by rectangular steps, with bounded integral and conditional anti-windup.
 * kp,kd,ki and limits are finite nonnegative values; dt>=0 seconds.
 * Outputs are body torque [N m], clipped componentwise. limited_axes optional.
 * Caller coordinates these bounds with actuator allocation. Atomic on failure. */
tinypid_status tinypid_step(tinypid_state *state,const double current[4],const double desired[4],
    const double rate_error[3],const double kp[3],const double kd[3],const double ki[3],
    const double feedforward[3],double dt,const double integral_limit[3],const double torque_limit[3],
    double torque[3],unsigned *limited_axes);
#ifdef __cplusplus
}
#endif
#endif
