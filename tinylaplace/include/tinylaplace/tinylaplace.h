#ifndef TINYLAPLACE_H
#define TINYLAPLACE_H
#include "tinypolynomial/tinypolynomial.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYLAPLACE_OK=0,TINYLAPLACE_INVALID_ARGUMENT,TINYLAPLACE_GEOMETRY,TINYLAPLACE_NO_SOLUTION,TINYLAPLACE_NO_CONVERGENCE,TINYLAPLACE_ARITHMETIC_ERROR} tinylaplace_status;
typedef struct {double position_m[3],velocity_m_s[3],range_m,range_rate_m_s,selected_radius_m;unsigned admissible_roots;} tinylaplace_result;
typedef struct {tinypolynomial_complex roots[8],scratch[16];} tinylaplace_workspace;
/* Instantaneous Laplace curvature method. LOS is a unit vector; its time
 * derivatives in s^-1/s^-2 must describe that same inertial direction. Observer
 * position/velocity/acceleration are SI and share the epoch/frame. mu m³/s².
 * Eighth-degree radius polynomial may have multiple physical roots: hint=0
 * selects largest positive-range root, positive hint selects nearest radius.
 * Outputs preserved on failure; scratch is caller-owned and separate. */
tinylaplace_status tinylaplace_estimate(double mu,const double direction[3],const double direction_rate[3],const double direction_acceleration[3],
    const double observer_position[3],const double observer_velocity[3],const double observer_acceleration[3],double radius_hint_m,
    tinylaplace_result *result,tinylaplace_workspace *workspace);
/* Three equally spaced directions at [-spacing,0,+spacing], normalized before
 * central differentiation. This is an approximate initial estimate; derivative
 * truncation error can dominate. Supply actual middle observer kinematics. */
tinylaplace_status tinylaplace_three(double mu,double spacing_s,const double directions[9],const double observer_position[3],
    const double observer_velocity[3],const double observer_acceleration[3],double radius_hint_m,tinylaplace_result *result,tinylaplace_workspace *workspace);
#ifdef __cplusplus
}
#endif
#endif
