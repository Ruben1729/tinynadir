#ifndef TINYPD_H
#define TINYPD_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYPD_OK=0,TINYPD_INVALID_ARGUMENT,TINYPD_DEGENERATE,TINYPD_ARITHMETIC_ERROR} tinypd_status;
/* Both q are scalar-first reference->body Hamilton quaternions. rate_error is
 * physical measured minus desired rate already expressed in current body axes.
 * tau=2*kp*sign(q_error.w)*q_error.xyz-kd*rate_error+feedforward,
 * q_error=q_current*inverse(q_desired). kp,kd are nonnegative per-axis gains;
 * optional feedforward is body torque, NULL means zero. Returns requested body
 * torque [N m]; actuator allocation/limits belong outside this controller. */
tinypd_status tinypd_torque(const double current[4],const double desired[4],const double rate_error[3],
    const double kp[3],const double kd[3],const double feedforward[3],double torque[3]);
#ifdef __cplusplus
}
#endif
#endif
