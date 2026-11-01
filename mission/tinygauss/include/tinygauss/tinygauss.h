#ifndef TINYGAUSS_H
#define TINYGAUSS_H
#include "tinypolynomial/tinypolynomial.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYGAUSS_OK=0,TINYGAUSS_INVALID_ARGUMENT,TINYGAUSS_GEOMETRY,TINYGAUSS_NO_SOLUTION,TINYGAUSS_NO_CONVERGENCE,TINYGAUSS_ARITHMETIC_ERROR} tinygauss_status;
typedef struct {double position_m[3],velocity_m_s[3],ranges_m[3],selected_radius_m;unsigned admissible_roots;} tinygauss_result;
typedef struct {tinypolynomial_complex roots[8],scratch_roots[16];double linear[12];} tinygauss_workspace;
/* Classical Gauss initial estimate: eighth-degree radius polynomial, truncated
 * f/g range coefficients, and Gibbs velocity. Not an exact fitted trajectory.
 * tau1<0<tau3 are seconds relative to middle epoch, |tau_i|<=1e6. LOS and
 * observer arrays contain three consecutive inertial xyz triples (directions
 * normalized internally; observer positions in metres). mu in m³/s².
 * hint=0 selects largest admissible real positive root; positive radius hint
 * selects nearest admissible radius. Every selected range must be positive.
 * Output is preserved on failure; caller workspace must be separate. */
tinygauss_status tinygauss_estimate(double mu,double tau1,double tau3,const double line_of_sight[9],
    const double observer_position_m[9],double radius_hint_m,tinygauss_result *result,tinygauss_workspace *workspace);
#ifdef __cplusplus
}
#endif
#endif
