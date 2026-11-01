#ifndef TINYCOMPLEMENTARY_H
#define TINYCOMPLEMENTARY_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYCOMPLEMENTARY_MAX_OBSERVATIONS 64u
typedef enum {TINYCOMPLEMENTARY_OK=0,TINYCOMPLEMENTARY_INVALID_ARGUMENT,TINYCOMPLEMENTARY_DEGENERATE,TINYCOMPLEMENTARY_ARITHMETIC_ERROR} tinycomplementary_status;
typedef struct {double q[4],bias[3];} tinycomplementary_state;
typedef struct {double proportional_gain,integral_gain,bias_limit;} tinycomplementary_parameters;
tinycomplementary_status tinycomplementary_initialize(tinycomplementary_state *state,const double q[4],const double bias[3]);
/* Mahony-style explicit observer for arbitrary reference/body unit directions.
 * q maps reference->body; gyro physical body rate [rad/s], dt>0 [s]. Arrays
 * hold count triples, count 0..64. NULL weights mean equal weights; nonnegative
 * weights are normalized. Zero count/all-zero weights is gyro-only prediction.
 * error=sum(weight*(measured cross predicted)), bias_dot=-ki*error,
 * corrected_rate=gyro-new_bias+kp*error. Uses normalized first-order quaternion
 * integration, not tinygyro's exact held-rate step. bias is component-limited.
 * No gravity/accelerometer direction is assumed. Failure preserves state. */
tinycomplementary_status tinycomplementary_step(tinycomplementary_state *state,const tinycomplementary_parameters *parameters,
    const double gyro[3],double dt,size_t count,const double *reference,const double *body,const double *weights);
#ifdef __cplusplus
}
#endif
#endif
