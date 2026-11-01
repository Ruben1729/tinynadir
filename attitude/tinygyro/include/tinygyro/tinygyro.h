#ifndef TINYGYRO_H
#define TINYGYRO_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYGYRO_OK=0,TINYGYRO_INVALID_ARGUMENT,TINYGYRO_DEGENERATE,TINYGYRO_ARITHMETIC_ERROR} tinygyro_status;
/* Scalar-first Hamilton q maps inertial/reference to body coordinates.
 * gyro is physical body angular rate relative to the reference, expressed in
 * body axes [rad/s]. Optional bias uses the same axes/units; NULL means zero.
 * qdot = -0.5*[0,gyro-bias] * q. step exactly integrates held corrected body
 * rate for dt>=0 seconds; it does not integrate arbitrary varying rates.
 * Failed operations preserve outputs; outputs may alias q. */
tinygyro_status tinygyro_derivative(const double q[4],const double gyro[3],const double bias[3],double qdot[4]);
tinygyro_status tinygyro_step(const double q[4],const double gyro[3],const double bias[3],double dt,double out[4]);
#ifdef __cplusplus
}
#endif
#endif
