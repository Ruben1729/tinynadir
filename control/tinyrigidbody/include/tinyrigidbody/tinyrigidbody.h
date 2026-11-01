#ifndef TINYRIGIDBODY_H
#define TINYRIGIDBODY_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYRIGIDBODY_OK=0,TINYRIGIDBODY_INVALID_ARGUMENT,TINYRIGIDBODY_INVALID_INERTIA,TINYRIGIDBODY_ARITHMETIC_ERROR} tinyrigidbody_status;
typedef struct {double factor[9],scratch[12];} tinyrigidbody_workspace;
/* Reference->body Hamilton q, body physical angular rate [rad/s], row-major
 * symmetric positive-definite inertia [kg m^2], applied body torque [N m],
 * optional body wheel momentum [N m s] (NULL zero). The torque already includes
 * reaction-wheel motor reaction; do not subtract its momentum derivative twice.
 * J*omega_dot=torque-omega cross (J*omega+wheel_momentum).
 * qdot follows tinygyro's physical sign. Outputs distinct; atomic on failure. */
tinyrigidbody_status tinyrigidbody_derivative(const double q[4],const double omega[3],const double inertia[9],
    const double applied_torque[3],const double wheel_momentum[3],double qdot[4],double omega_dot[3],tinyrigidbody_workspace *workspace);
#ifdef __cplusplus
}
#endif
#endif
