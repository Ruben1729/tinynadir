#ifndef TINYROTATION_H
#define TINYROTATION_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    TINYROTATION_OK = 0,
    TINYROTATION_INVALID_ARGUMENT,
    TINYROTATION_DEGENERATE,
    TINYROTATION_NOT_A_ROTATION,
    TINYROTATION_ARITHMETIC_ERROR
} tinyrotation_status;

/* Scalar-first Hamilton q=[w,x,y,z]. C(q) maps reference vector coordinates to
 * body coordinates; C(p*q)=C(p)*C(q). Matrices are row-major. Input quaternions
 * and axes are normalized robustly. Arrays may alias except separate outputs
 * of to_axis_angle. Failed calls preserve outputs. Angles are radians. */
tinyrotation_status tinyrotation_normalize(const double q[4], double out[4]);
tinyrotation_status tinyrotation_multiply(const double p[4], const double q[4], double out[4]);
tinyrotation_status tinyrotation_inverse(const double q[4], double out[4]);
tinyrotation_status tinyrotation_to_dcm(const double q[4], double dcm[9]);
tinyrotation_status tinyrotation_from_dcm(const double dcm[9], double q[4]);
tinyrotation_status tinyrotation_apply(const double q[4], const double reference[3], double body[3]);
tinyrotation_status tinyrotation_from_axis_angle(const double axis[3], double angle, double q[4]);
tinyrotation_status tinyrotation_to_axis_angle(const double q[4], double axis[3], double *angle);
tinyrotation_status tinyrotation_from_rotvec(const double vector[3], double q[4]);
tinyrotation_status tinyrotation_to_rotvec(const double q[4], double vector[3]);
/* Intrinsic ZYX / extrinsic XYZ: C=Rz(yaw)*Ry(pitch)*Rx(roll).
 * to_euler_zyx returns [yaw,pitch,roll], pitch in [-pi/2,pi/2]. At gimbal
 * lock it selects roll=0; individual yaw/roll are not observable there. */
tinyrotation_status tinyrotation_from_euler_zyx(double yaw, double pitch, double roll, double q[4]);
tinyrotation_status tinyrotation_to_euler_zyx(const double q[4], double angles[3]);
/* Shortest-arc interpolation for fraction in [0,1]. */
tinyrotation_status tinyrotation_slerp(const double from[4], const double to[4], double fraction, double out[4]);

#ifdef __cplusplus
}
#endif
#endif
