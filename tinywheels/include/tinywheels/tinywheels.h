#ifndef TINYWHEELS_H
#define TINYWHEELS_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYWHEELS_MAX_WHEELS 8u
typedef enum {TINYWHEELS_OK=0,TINYWHEELS_INVALID_ARGUMENT,TINYWHEELS_NO_AUTHORITY,TINYWHEELS_ARITHMETIC_ERROR} tinywheels_status;
typedef struct {double gram[9],values[3],vectors[9],scratch[18];} tinywheels_workspace;
/* count 1..8, axes are count consecutive unit body-axis triples. wheel torque
 * is applied to rotors; body reaction is -sum(axis_i*u_i). Compute minimum-norm
 * unconstrained allocation with a symmetric-eigen pseudoinverse, then clip
 * torque/speed authority; no iterative redistribution after clipping occurs.
 * enabled may be NULL; speed [rad/s], limits nonnegative [N m,rad/s].
 * Residual=requested_body_torque-achieved_body_torque exposes lost authority.
 * Output arrays must be distinct, workspace separate; failure preserves outputs. */
tinywheels_status tinywheels_allocate(size_t count,const double *axes,const double requested[3],
    const double *speed,const double *torque_limit,const double *speed_limit,const unsigned char *enabled,
    double *wheel_torque,double residual[3],tinywheels_workspace *workspace);
#ifdef __cplusplus
}
#endif
#endif
